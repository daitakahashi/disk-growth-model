
// [[Rcpp::plugins(cpp11)]]
#include <Rcpp.h>

#include <algorithm>
#include <vector>
#include <tuple>

#include <iostream>

struct Region {
  using CellType = char;  
  using Matrix = std::vector<CellType>;
  static const CellType exist  = 1;
  static const CellType absent = 0;
};

const Region::CellType Region::exist;
const Region::CellType Region::absent;

struct IntensityDistribution;
struct LocalPopulation;
struct RegionalPopulation;

struct IntensityDistribution {
  using Matrix = Rcpp::NumericMatrix;
  Matrix base;
  Matrix dist;
  double total;
  const int ntotal;
  int nactive;
  
  int total_cell_count() const {
    return ntotal;
  }
  int active_cell_count() const {
    return nactive;
  }

  IntensityDistribution(Matrix& dist) :
    base(dist), dist(dist), total(Rcpp::sum(dist)),
    ntotal(dist.size()),
    nactive(ntotal) {}

  int exclude_covered_area(const RegionalPopulation& covered_area);
  void update_total() {
    total = Rcpp::sum(dist);
  }
  void reset() {
    dist = base;
    update_total();
    nactive = dist.size();
  }
};

struct LocalPopulation {
  const double x;
  const double y;
  double rad;
  int age;
  bool active;
  
  LocalPopulation(const double x,
                  const double y,
                  const double r,
		  const int age)
    : x(x), y(y), rad(r), age(age), active(true) {}
  LocalPopulation(const double x,
                  const double y,
                  const double r)
    : x(x), y(y), rad(r), age(0), active(true) {}
};

struct Statistics {
  using nameType = std::vector<std::string>;
  using dataType = std::vector<double>;
  nameType names;
  dataType data;
  
  void assign(const std::string key, const double val) {
    names.push_back(key);
    data.push_back(val);
  }
};

bool name_exists(const Rcpp::CharacterVector& v, const char* name) {
  return std::find(v.begin(), v.end(), name) != v.end();
}

struct RegionalPopulation {
  using Matrix = Region::Matrix;
  using LocalPopulations = std::vector<LocalPopulation>;
  
  const int spatial_resolution;
  Matrix region;
  const size_t total_cells;
  LocalPopulations pop;
  size_t total_covered_cells = 0;
  size_t n_active_populations = 0;

  RegionalPopulation(const Rcpp::NumericMatrix initial_population,
                     const int resolution,
		     const double scale = 1)
    : spatial_resolution(resolution),
      region(resolution*resolution),
      total_cells(region.size())
  {
    const int n_population = initial_population.nrow();
    pop.reserve(n_population);

    std::fill(region.begin(), region.end(), Region::absent);

    Rcpp::NumericMatrix::ConstColumn X = initial_population(Rcpp::_, 0);
    Rcpp::NumericMatrix::ConstColumn Y = initial_population(Rcpp::_, 1);
    Rcpp::NumericMatrix::ConstColumn radius = initial_population(Rcpp::_, 2);
    
    if (initial_population.ncol() == 4) {
      // initial-population matrix also specifies ages of local populations
      Rcpp::NumericMatrix::ConstColumn ages
	= initial_population(Rcpp::_, 3);
      
      for (int i = 0; i < n_population; i++) {
	// convert from 1-base cell-centor index
	pop.emplace_back((X[i] - 0.5)*scale, (Y[i] - 0.5)*scale,
			 radius[i]*scale,
			 ages[i]);
      }
    } else {
      for (int i = 0; i < n_population; i++) {
	// convert from 1-base cell-centor index
	pop.emplace_back((X[i] - 0.5)*scale, (Y[i] - 0.5)*scale,
			 radius[i]*scale);
      }
    }

    n_active_populations = n_population;
    // draw();
  }
  
  void draw();
  void increment_step(const double g, const double s,
                      const IntensityDistribution& pickupIntensities,
                      const IntensityDistribution& releaseProbabilities);
  
  void deactivate_inner_populations();
  void deactivate_all();

  double total_covered_proportion () const {
    return (double)total_covered_cells/total_cells;
  }

  void register_statistics(Statistics& stat);
  
  size_t mark_inactive(const double margin);
  
  void growth(const double g, const double s);
  void birth(const IntensityDistribution& pickupIntensities,
             const IntensityDistribution& releaseProbabilities);

  Rcpp::IntegerVector
  population_on_region(Rcpp::IntegerMatrix& mask);

  void remove_populations(Rcpp::IntegerVector& mask_removal);
};

struct SimulationSettings {
  const size_t n_iterations;
  const bool stop_when_complete;
  const bool print_progress;
  const bool include_populations_in_result;
  const double target_proportion;
  const bool remove_detected_colonies;
  const double minimum_detectable_size;

  SimulationSettings(const size_t _n_iterations,
		     const bool _stop_when_complete,
		     const bool _print_progress,
		     const bool _include_populations_in_result,
		     const double _target_proportion = 1,
                     const bool _remove_detected_colonies = false,
                     const double _minimum_detectable_size = 0)
    : n_iterations(_n_iterations),
      stop_when_complete(_stop_when_complete),
      print_progress(_print_progress),
      include_populations_in_result(_include_populations_in_result),
      target_proportion(_target_proportion),
      remove_detected_colonies(_remove_detected_colonies),
      minimum_detectable_size(_minimum_detectable_size) {}
};

struct Result {
  RegionalPopulation::LocalPopulations pop;
  Statistics stat;
};


template <typename T>
static inline
T square(const T x) {
  return x*x;
}

void RegionalPopulation::draw() {
  const int n = spatial_resolution;
  const LocalPopulations& p = pop;
  const int nccl = p.size();
  auto region_p = region.begin();

  for (int k = 0; k < nccl; k++) {
    const double x = p[k].x;
    const double y = p[k].y;
    const double r = p[k].rad;
    const bool area_active = p[k].active;

    if (r <= 0) continue;
    
    if (area_active) {
      const double r2 = square(r);
      // compute a bbox of a k-th disk
      const int minx = std::max(0, (int)std::floor(x - r));
      const int miny = std::max(0, (int)std::floor(y - r));
      const int endx = std::min(n, (int)std::ceil(x + r) + 1);
      const int endy = std::min(n, (int)std::ceil(y + r) + 1);
      
      // iterate over the bbox
      for (int j = miny; j < endy; j++) {
        const double discriminant = r2 - square(y - j);
        if (discriminant >= 0) {
          const double d = sqrt(discriminant);
          const int init_i = std::min(n - 1,
                                      std::max(minx, (int)std::floor(x - d)));
          const int end_i  = std::max(init_i,
                                      std::min(endx, (int)std::ceil(x + d)));
          
          std::fill(region_p + j*n + init_i,
                    region_p + j*n + end_i,
                    Region::exist);
        }
      }
    }
  }

  total_covered_cells
    = std::accumulate(region.begin(), region.end(), 0);
  
  return;
}

// *** ASSUME "coveredArea" IS REUSED (WITHOUT RESET) AT EVERY TIME STEP. ***
// Then we can skip drawing of an already drawn population if surrounding area
// of the population has also occupied (i.e., no room for the lateral growth).
static
size_t update_active_marks(Rcpp::IntegerVector& active_mark,
                           const RegionalPopulation::LocalPopulations& pop,
                           const RegionalPopulation::Matrix& region,
                           const size_t spatial_resolution,
                           const double margin = 1.5) {  
  const size_t n = spatial_resolution;
  const size_t n_total_population = pop.size();
  auto region_p = region.begin();
  size_t n_deactivated = 0;

  // A local population that is completelly surrounded by other
  // local populations is to be deactivated, since following expansions
  // of the population no longer affect total area size
  // covered by local populations.
  for (size_t k = 0; k < n_total_population; k++) {
    const double x = pop[k].x;
    const double y = pop[k].y;
    const double r = pop[k].rad + margin;
    
    if (r <= 0 && active_mark[k] == 1) {
      if (region[floor(y)*n + floor(x)] == Region::exist) {
        active_mark[k] = 0;
        n_deactivated++;
      }
    } else if (active_mark[k] == 1) {
      const double r2 = square(r);
      const int minx = std::max(0, (int)std::floor(x - r));
      const int endx = std::min((int)n, (int)std::ceil(x + r) + 1);
      const int miny = std::max(0, (int)std::floor(y - r));
      const int endy = std::min((int)n, (int)std::ceil(y + r) + 1);
      
      // iterate over the bbox
      for (size_t j = miny; j < (size_t)endy; j++) {
        const double discriminant = r2 - square(y - j);
        
        if (discriminant >= 0) {
          const double d = sqrt(discriminant);
          const int init_i = std::min((int)n - 1,
                                         std::max(minx, (int)std::floor(x - d)));
          const int end_i  = std::max(init_i,
                                         std::min(endx, (int)std::ceil(x + d)));
          auto iter_first = region_p + j*n + init_i;
          auto iter_last = region_p + j*n + end_i;
          
          if (std::find(iter_first, iter_last, Region::absent) != iter_last) {
            // found an open cell; keep the current population to be active
            active_mark[k] = 1;
            goto next_population;
          }
        }
      }
      // did not find any open cells around this population.
      active_mark[k] = 0;
      n_deactivated++;
    }
    next_population:
    {/* do nothing more on this population and continue to the next one */}
  }
  return n_deactivated;
}


size_t RegionalPopulation::mark_inactive(const double margin = 1.5) {
  const size_t n = pop.size();
  Rcpp::IntegerVector active_mark(n);
  for (size_t i = 0; i < n; i++) {
    active_mark[i] = pop[i].active;
  }

  const size_t n_deactivated = update_active_marks(active_mark,
                                                   pop,
                                                   region,
                                                   spatial_resolution,
                                                   margin);
  
  for (size_t i = 0; i < n; i++) {
    pop[i].active = active_mark[i];
  }
  n_active_populations = Rcpp::sum(active_mark);

  return n_deactivated;
}

void RegionalPopulation::growth(const double g, const double s) {
  for (auto& x : pop) {
    x.age++;
    x.rad += R::rnorm(g, g*s);
    if (x.rad < 0) {
      x.rad = 0;
    }
  }
  return;
}

static
double masked_sum(const Region::Matrix& mask,
                  const IntensityDistribution::Matrix& v) {
  const std::size_t size = mask.size();
  double result = 0;

  for (std::size_t i = 0; i < size; i++) {
    if (mask[i] == Region::exist) {
      result += v[i];
    }
  }
    
  return result;
}

void RegionalPopulation::birth(const IntensityDistribution& pickup_intensities,
                               const IntensityDistribution& destination_probs) {
  const double average_num_disperser
    = (masked_sum(region, pickup_intensities.dist)
       *destination_probs.total);
  const int n_dispersers = R::rpois(average_num_disperser);
  
  if (n_dispersers > 0) {
    const Rcpp::NumericMatrix& dist = destination_probs.dist;
    const Rcpp::IntegerVector ixs = Rcpp::sample(dist.size(),
                                                 n_dispersers,
                                                 true,
                                                 dist,
                                                 false);
    
    const int nr = spatial_resolution;
    const int dp_nr = dist.nrow();
    const double dp_fac = (double)nr/dp_nr;

    for (int i = 0; i < n_dispersers; i++) {
      const int ix = ixs[i];
      const double x = dp_fac*((ix % dp_nr) + R::runif(0, 1));
      const double y = dp_fac*((ix/dp_nr)   + R::runif(0, 1));
      
      const int x_gp = std::max(std::min(nr - 1, (int)floor(x)), 0);
      const int y_gp = std::max(std::min(nr - 1, (int)floor(y)), 0);
      
      if (region[y_gp*spatial_resolution + x_gp] < 1) {
        pop.emplace_back(x, y, 0.0);
	n_active_populations++;
      }
    }
  }
  
  return;
}

void RegionalPopulation::increment_step(const double g,
                                        const double s,
                                        const IntensityDistribution& pickup_intensities,
                                        const IntensityDistribution& release_probabilities) {
  
  this->growth(g, s);
  this->draw();
  this->birth(pickup_intensities, release_probabilities);

  return;
}

void RegionalPopulation::register_statistics(Statistics& stat) {
  stat.assign("noccupied", total_covered_cells);
  stat.assign("npopulations", pop.size());
  size_t n_new_populations = 0;
  for (auto& x : pop) {
    if (x.age == 0) {
      n_new_populations++;
    }
  }
  stat.assign("npopulations.new", n_new_populations);
}

std::ostream& operator<<(std::ostream& ostr,
			 const RegionalPopulation& rpop) {
  ostr << "npop: " << rpop.pop.size() << " "
       << "(" << rpop.n_active_populations << " active), "
       << "covered: " << rpop.total_covered_cells << " / " << rpop.region.size();
  return ostr;
}


void RegionalPopulation::deactivate_inner_populations() {
  
  this->mark_inactive();

  return;
}

void RegionalPopulation::deactivate_all() {
  
  for (auto& x : pop) {
    x.active = false;
  }
  
  n_active_populations = 0;
  return;
}

int IntensityDistribution::exclude_covered_area(const RegionalPopulation& rpop) {
  const Region::Matrix& region = rpop.region;
  const int n_region = rpop.spatial_resolution;
  const int n_dist   = dist.nrow();
  const double ratio = (double)n_region/n_dist;
  int n_deactivate   = 0;

  for (int j = 0; j < n_dist; j++) {
    for (int i = 0; i < n_dist; i++) {
      if (dist(i, j) > 0) {
        const int rj_end = std::min(n_region, (int)std::ceil((j+1)*ratio));
        for (int rj = std::floor(j*ratio); rj < rj_end; rj++) {
          const int x_begin = std::floor(i*ratio);
          const int x_end   = std::min(n_region, (int)std::ceil((i+1)*ratio));
          auto iter_begin   = region.begin() + rj*n_region + x_begin;
          auto iter_end     = region.begin() + rj*n_region + x_end;
          
          if (std::find(iter_begin, iter_end, Region::absent) != iter_end) {
            // at least a cell in the block is open
            // keep (i,j) to be a dispersal destination
            goto next_cell;
          }
        }
        // all cells in the block is covered by populations.
        // the establishment of a disperser comming to an interval
        // (i, j) + ([0, 1), [0, 1))
        // always fail because everywhere in the interval is already occupied
        // therefore, exclude (i, j) from dispersal destinations
        dist(i, j) = 0;
        n_deactivate++;
      }
    next_cell:
      {/* do nothing more on this cell and continue to the next one */}
    }
  }

  this->update_total();
  nactive -= n_deactivate;
  
  return n_deactivate;
}

static
int deactivate_redundant_elements(RegionalPopulation& rpop,
                                  IntensityDistribution& disperser_release_prob,
                                  const SimulationSettings& settings) {
  
  rpop.deactivate_inner_populations();
  
  const int n_deactivated
    = disperser_release_prob.exclude_covered_area(rpop);
  
  if (settings.print_progress) {
    Rcpp::Rcout << "deactivate " << n_deactivated << " dispersal-target cells"
                << " (" << disperser_release_prob.active_cell_count()
                << " of " << disperser_release_prob.total_cell_count()
                << " cells currently active)"
                << std::endl;
  }

  return n_deactivated;
}

static
void print_population(const size_t i,
		      const RegionalPopulation& rpop,
		      const SimulationSettings& settings) {
  if (settings.print_progress) {
    Rcpp::Rcout << "t: " << i << ", "
		<< rpop
		<< std::endl;
  }
  return;
}

static
bool area_reaches_target_proportion(const RegionalPopulation& rpop,
				    const double target_proportion) {
  if (target_proportion == 1) {
    return rpop.total_covered_cells >= rpop.region.size();
  } else {
    return (double)rpop.total_covered_cells/rpop.region.size() >= target_proportion;
  } 
}

static
Rcpp::IntegerMatrix bwlabel(Rcpp::IntegerMatrix m) {
  Rcpp::Environment ebimage("package:EBImage");
  Rcpp::Function bwlabel = ebimage["bwlabel"];
  Rcpp::IntegerMatrix result = bwlabel(Rcpp::as<SEXP>(m));
  return result;
}

static
std::tuple<Rcpp::IntegerMatrix, size_t>
make_mask(Rcpp::IntegerMatrix labels,
                              const double threshold_area) {
  const int nrow = labels.nrow();
  const int ncol = labels.ncol();
  Rcpp::IntegerMatrix result(nrow, ncol);
  const int label_max = Rcpp::max(labels);
  size_t n_removed = 0;
  
  if (label_max > 0) {
    std::vector<int> sizes(label_max);
    std::fill(sizes.begin(), sizes.end(), 0);
    for (int x : labels) {
      if (x > 0) {
        sizes[x - 1]++;
      }
    }
  
    for (int i = 0; i < nrow*ncol; i++) {
      const int val = labels[i];
      if (val > 0 && sizes[val - 1] > threshold_area) {
        result[i] = 1;
      }
    }
    n_removed = std::count_if(sizes.begin(), sizes.end(),
			      [threshold_area](int x){return x > threshold_area;});
  }
  return std::make_tuple(result, n_removed);
}

Rcpp::IntegerVector
RegionalPopulation::population_on_region(Rcpp::IntegerMatrix& mask) {
  const size_t n = pop.size();
  Rcpp::IntegerVector result(n);
  std::fill(result.begin(), result.end(), 1);
  Matrix mask_mat(mask.size());
  std::copy(mask.begin(), mask.end(), mask_mat.begin());

  update_active_marks(result,
                      pop,
                      mask_mat,
                      spatial_resolution,
                      0);

  // active : the population locates outside of the masked region
  //  --> "1 - result" specifies populations inside of the region
  return 1 - result; 
}

void RegionalPopulation::remove_populations(Rcpp::IntegerVector& mask_removal) {
  const size_t current_n = pop.size();
  LocalPopulations new_population;
  
  for (size_t i = 0; i < current_n; i++) {
    if (mask_removal[i] != 1) {
      new_population.emplace_back(std::move(pop[i]));
    }
  }
  pop.clear();
  pop = std::move(new_population);
  
  const size_t new_n = pop.size();

  for (size_t i = 0; i < new_n; i++) {
    pop[i].active = true;
  }
  
  n_active_populations = pop.size();
  std::fill(region.begin(), region.end(), 
            Region::absent);
  
  return;
}

static
bool finish_step(RegionalPopulation& rpop,
                 std::vector<Result>& result_vec,
                 IntensityDistribution& disperser_release_prob,
                 const size_t i,
                 const SimulationSettings& settings) {
  Statistics stat;
  const RegionalPopulation::LocalPopulations
    current_populations(rpop.pop);
  
  rpop.register_statistics(stat);

  Rcpp::IntegerMatrix m(rpop.spatial_resolution,
                        rpop.spatial_resolution);
  std::copy(rpop.region.begin(), rpop.region.end(),
            m.begin());

  Rcpp::IntegerMatrix labels = bwlabel(m);
  stat.assign("nclusters", Rcpp::max(labels));

  if (settings.remove_detected_colonies) {
    std::tuple<Rcpp::IntegerMatrix, int> removal_plan
      = make_mask(labels,
		  settings.minimum_detectable_size);
    Rcpp::IntegerVector removed_populations
      = rpop.population_on_region(std::get<0>(removal_plan));

    stat.assign("nclusters.removed", std::get<1>(removal_plan));
    stat.assign("npopulations.removed", Rcpp::sum(removed_populations));
    stat.assign("area.removed", Rcpp::sum(std::get<0>(removal_plan)));

    rpop.remove_populations(removed_populations);

  } else {
    // try to deactivate local populations and dispersal targets once in 10 steps
    if ((i % 10) == 0) {
      deactivate_redundant_elements(rpop, disperser_release_prob, settings);
    }
  }

  result_vec.push_back(Result({current_populations, stat}));
  
  // all cells are covered. skip following iterations
  if (area_reaches_target_proportion(rpop, settings.target_proportion)
      && settings.stop_when_complete) {
    
    rpop.deactivate_all();
    
    if (settings.print_progress) {
      Rcpp::Rcout << "t: " << i
                  << ": the region is completely covered (skip followings)"
                  << std::endl;
    }
    return false;
  }

  // all populations are died
  if (rpop.pop.size() == 0) {
    return false;
  }
  return true;
}

static
Rcpp::NumericMatrix
convert_pop_to_matrix(const RegionalPopulation::LocalPopulations& pop) {
  const size_t n_populations = pop.size();
  Rcpp::NumericMatrix result(n_populations, 5);

  for (size_t i = 0; i < n_populations; i++) {
    result(i, 0) = pop[i].x + 0.5;
    result(i, 1) = pop[i].y + 0.5;
    result(i, 2) = pop[i].rad;
    result(i, 3) = pop[i].age;
    result(i, 4) = pop[i].active;
  }
  colnames(result) = Rcpp::CharacterVector::create("x", "y", "radius",
                                                   "age", "active");
  
  return result;
}

static
Rcpp::NumericMatrix convert_stats_to_matrix(const std::vector<Result>& v) {
  const size_t len = v.size();
  const size_t n_stats = v[0].stat.data.size();
  Rcpp::NumericMatrix stat_matrix(len, n_stats);

  for (size_t i = 0; i < len; i++) {
    for (size_t j = 0; j < n_stats; j++) {
      stat_matrix(i, j) = v[i].stat.data[j];
    }
  }
  Rcpp::CharacterVector names(n_stats);
  std::copy(v[0].stat.names.begin(), v[0].stat.names.end(),
            names.begin());
  colnames(stat_matrix) = names;
  return stat_matrix;
}

static
Rcpp::List convert_output_to_list(const std::vector<Result>& v,
                                  const SimulationSettings& settings) {
  const size_t len = v.size();
  
  if (len == 0) {
    return Rcpp::List::create();
  }

  Rcpp::NumericMatrix stat_matrix = convert_stats_to_matrix(v);
  
  if (settings.include_populations_in_result) {
    Rcpp::List result_pops(len);
    for (size_t i = 0; i < len; i++) {
      result_pops[i] = convert_pop_to_matrix(v[i].pop);
    }
    return Rcpp::List::create(Rcpp::Named("populations") = result_pops,
                              Rcpp::Named("stats") = stat_matrix);
  } else {
    return Rcpp::List::create(Rcpp::Named("stats") = stat_matrix);
  }
}

// [[Rcpp::export]]
Rcpp::List runSimulation(const int nIterations,
                         const double g,
                         const double s,
                         Rcpp::NumericMatrix disperserPickupIntensities,
                         Rcpp::NumericMatrix disperserReleaseProbabilities,
                         const Rcpp::NumericMatrix initialPopulation,
                         const bool stopWhenComplete = false,
                         const bool print_progress = true,
                         const bool include_populations_in_result = true,
                         const double target_proportion = 1,
                         const bool remove_detected_colonies = false,
                         const double minimum_detectable_size = 0) {

  const SimulationSettings settings(nIterations,
				    stopWhenComplete,
				    print_progress,
				    include_populations_in_result,
				    target_proportion,
                                    remove_detected_colonies,
                                    minimum_detectable_size);

  const int rasterizeResolution = disperserPickupIntensities.nrow();
  IntensityDistribution disperser_pickup_intensity(disperserPickupIntensities);

  Rcpp::NumericMatrix disperser_release_prob_matrix
    = disperserReleaseProbabilities/Rcpp::sum(disperserReleaseProbabilities);
  IntensityDistribution disperser_release_prob(disperser_release_prob_matrix);

  std::vector<Result> result_vec;
  
  // instanciate an initial state (i == 0) and export it to a result vector
  RegionalPopulation rpop(initialPopulation, rasterizeResolution);

  // deactivate_redundant_elements(rpop, disperser_release_prob, print_progress);
  rpop.draw();

  finish_step(rpop, result_vec, disperser_release_prob, 0, settings);

  print_population(0, rpop, settings);

  // main loop (i starts from 1)
  for (int i = 1; i < nIterations + 1; i++) {
    // check for user interrupt every 10 steps
    if (i % 10 == 0) {
        Rcpp::checkUserInterrupt();
    }    

    // main routine
    rpop.increment_step(g, s, disperser_pickup_intensity, disperser_release_prob);
    
    // print model status
    print_population(i, rpop, settings);

    const bool continue_next = finish_step(rpop, result_vec,
                                           disperser_release_prob,
                                           i,
                                           settings);
    if (!continue_next) {
      break;
    }
  }
  
  return Rcpp::List::create(convert_output_to_list(result_vec, settings),
                            convert_pop_to_matrix(rpop.pop));
}

// [[Rcpp::export]]
Rcpp::NumericMatrix drawOnMatrixC(const int n,
                                  Rcpp::NumericMatrix ccl,
				  const double scale) {
  RegionalPopulation rpop(ccl, n, /* scale = */ scale);
  Rcpp::NumericMatrix result(n, n);

  rpop.draw();
  // in general, rpop.region is not an instance of NumericMatrix.
  // we need to convert the rpop.region explicitly
  std::copy(rpop.region.begin(), rpop.region.end(),
            result.begin());
  return result;
}

// [[Rcpp::export]]
Rcpp::IntegerVector populationOnRegion(const int n,
				       Rcpp::NumericMatrix ccl,
				       Rcpp::NumericMatrix image) {
  const int m = image.nrow();
  const int npop = ccl.nrow();
  RegionalPopulation rpop(ccl, m, /* scale = */ (double)m/n);
  Rcpp::IntegerVector result(npop);

  std::copy(image.begin(), image.end(), rpop.region.begin());
  for (auto& x : rpop.region) {
    x = 1 - x;
  }
  
  rpop.mark_inactive(/* margin = */ 0.0);

  for (int i = 0; i < npop; i++) {
    result[i] = rpop.pop[i].active;
  }
  
  return result;
}
