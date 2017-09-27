# --------------------------------------------------------------------------- #
# nCOP - Network Coverage of Patients
# --------------------------------------------------------------------------- #
#
# Author: Borislav Hristov (borislav@cs.princeton.edu)
#
# I. Input 
#
# There are two required inputs:
# 1) a network file and 2) a mutational file containing genes with list of
# individuals that have variants in these genes.
# Additionally, the user may provide:
# 3) a weight file specifying a weight for each node in the network
# 4) a value for alpha (in which case the program skips the step selecting 
# alpha and uses the user specified value)
# 5) output prefix which is used in the beginning of the name of 
# the output file
#
# II. Output
#
# output_prefix_results.txt is written in the Outputs directory. The file
# contain a list of candidate genes ranked by how frequently they appear
# in the randomized runs.
#
# III. How to run
#
# 1. To run with basic inputs:
# ./run_nCOP network_file.txt mutational_file.txt 
#
# 2. If you want to specify any additional parameter add "param=value"
# ./run_nCOP network_file.txt mutational_file.txt weights=weights_file.txt alpha=0.5 output_prefix=my_output
#
# Note that nCOP is implemented in Ruby. You need an installed Ruby which can
# be simply done via: sudo apt-get install ruby-full
#
# IV. Input File Formats
#
# 1. Network file: each line specifies an edge, white space delimited:
# GENE_ID GENE_ID
#
# 2. Mutational file: each line is white space delimited
# GENE_ID INDIVIDUAL_1 INDIVIDUAL_3
# GENE_ID INDIVIDUAL_5
#
# 3. Weights file: 
# GENE_ID WEIGHT
# i.e TP53 0.59
#
# --------------------------------------------------------------------------- #


# --------------------------------------------------------------------------- #
# Read User Input
# --------------------------------------------------------------------------- #
def read_user_input
  if ARGV.size < 2
    puts "Insufficient inputs provided."
    exit
  
  elsif ARGV.size == 2
    $network       = ARGV[0]
    $dataMut       = ARGV[1]
    $weights_file  = "none"
    $alpha         = nil
    $cancer        = "nCOP_out"
  
  elsif ARGV.size > 2
    $network       = ARGV[0]
    $dataMut       = ARGV[1]
    $weights_file  = "none"
    $alpha         = nil
    $cancer        = "nCOP_out"
    
    (2..(ARGV.size - 1)).each do |i|
      index = ARGV[i].index("=")

      if index.nil?
        puts "Illegal input format. Please use param=value without white spaces if adding any non-required parameters, i.e alpha=0.4"
        exit
      end
      
      param = ARGV[i][0..(index - 1)]
      value = ARGV[i][(index+1)..-1]
      
      case param
      when "alpha"
        $alpha = value.to_f
      when "weights"
        $weights_file = value
      when "output_prefix"
        $cancer = value
      else
        puts "Wrong parameter specified. Please use one of the alpha=, weights=, or output_prefix=."
        exit
      end
    end
  end
end

def run_according_to_user_choices
  # read the network and the mutation input file
  read_network
  read_data_mutations
  
  # did the user specify a weights file?
  if $weights_file == "none"
    assign_default_weights
  else
    read_weights
  end
  
  # compute the normalization constant
  compute_normalization_cost
  
  # did the user spcify alpha?
  if $alpha.nil?
    run_range_alpha
    pick_alpha
  end
  
  # run nCOP
  run_nCOP_given_alpha
end

# --------------------------------------------------------------------------- #
# Read Input Files
# --------------------------------------------------------------------------- #
def read_network
  $g = RubGraph.new
  File.foreach($network) do |line|
    words = line.delete("\n").split()
    $g.add_edge(words[0], words[1])
  end
  $in_network = Set.new($g.nodes.keys)
  
  puts "Read network file: #{$network}"
  puts "Number of nodes = #{$g.nodes.keys.size}"
  puts "Number of edges = #{$g.num_edges}"
  puts ""  
end

def read_data_mutations
  $all_genes = Hash.new {|hash, key| hash[key] = Set.new() }
  $all_pat = Hash.new {|hash, key| hash[key] = Set.new() }
  genes_not_in_network = Set.new()
  
  # read the file
  File.foreach($dataMut) do |line|
    words = line.delete("\n").split()
    gene = words[0]

    # skip genes not present in the network
    unless $g.nodes.has_key?(gene)
      genes_not_in_network.add(gene)
      next
    end
    
    # assemble the hashes
    words.each_with_index do |pat, index|
      next if index == 0
      $all_genes[gene] << pat
      $all_pat[pat] << gene
    end
  end

  puts "Read mutation file: #{$dataMut}"
  puts "Genes not found in the network and excluded = #{genes_not_in_network.size}"
  puts "Number of mutated genes present = #{$all_genes.size}"
  puts "Number of patients = #{$all_pat.size}"
  puts ""
end

def read_weights
  $all_genes_weight = {}
  File.foreach($weights_file) do |line|
    words = line.delete("\n").split()
    $all_genes_weight[words[0]] = words[1].to_f
  end
  
  puts "Read weights file: #{$weights_file}"
  puts ""
end

def assign_default_weights
  $all_genes_weight = {}
  $g.nodes.keys.each { |gene| $all_genes_weight[gene] = 1 }
  
  puts "No weights file specified. Using default value: each node has a weight of 1."
  puts ""
end

# --------------------------------------------------------------------------- #
# Compute Normalization Cost
# --------------------------------------------------------------------------- #

def compute_normalization_cost
  puts "Computing normalization constant..."
  start_nrm = Time.now
  
  # set alpha and beta to the extremes
  $ALPHA = 1; $BETA = 0
  
  # keep track of normalization costs and fraction of covered patients at each iteration
  hsh_iter, hsh_nrmlz_costs = {}, {}
  
  # Iterate 10 times
  10.times do |iter|
    $iter = iter
    
    $all_avail_genes = deep_copy_hash $all_genes
    $all_avail_pat   = deep_copy_hash $all_pat
    $test_pat = Hash.new {|hash, key| hash[key] = Set.new() }
    pick_patients
  
    # run nCOP heuristic
    cover = exec_search "nCOP_heruistic"
  
    # compute the total weight of the graph
    total_normlz = 0
    cover.each { |e| total_normlz += $all_genes_weight[e] }
    
    # store the fraction of covered patients and the corresponding total weight
    hsh_iter[iter] = (get_num_covered($avail_pat,cover).to_f/$avail_pat.keys.size).round(2)
    hsh_nrmlz_costs[iter] = total_normlz.round(4)
  end
  
  # the maximum fraction of covered patients is
  max_frac = hsh_iter[sort_hash(hsh_iter)[0]]  
  
  # average over these normalization weights for which the cover is at least 90% from the max_frac
  arr_nrmlz = []
  hsh_nrmlz_costs.each_pair { |i,tn| arr_nrmlz << tn if hsh_iter[i] > 0.9*max_frac}
  $normalization_const = arr_nrmlz.avg.round(2)
  
  #puts "Computed normalization const = #{$normalization_const}"
  #puts "time computing normalization const = #{proper_time_since(start_nrm)}"
  puts "Done."
  puts ""
end

# --------------------------------------------------------------------------- #
# Starting Gene Selection
# --------------------------------------------------------------------------- #

# get the starting gene
def get_start_gene
  h, sorted_h, top = {}, {}, 0 
  
  if $weights_file == "none"
    # construct a hash using the number of patients each node covers
    $genes.each { |gene, pat| h[gene] = (10.to_f/$all_genes[gene].size).round(5) }
  else    
    # construct a hash using the weight of the nodes
    $genes.each { |gene, pat| h[gene] = $all_genes_weight[gene] }
  end
  
  # sort the given hash and select only its top genes
  h.values.sort.each do |v|
    h.map{ |k,v1| (v1==v) ? k : nil }.compact.each { |e| sorted_h[e] = (1/v).round(5) }
    top += 1
    break if top >= $top_limit
  end   
  
  # pick probabilistically a starting gene based on the prob specified by the hash

  # initialize a pick up
  #require 'pickup'
  pickup = Pickup.new(sorted_h)
    
  # make sure we don't pick bad starting genes
  num_try = 0
  while true do
    # we pick this gene
    start_gene = pickup.pick
    $picked[start_gene] += 1
    # repeat if it is a bad one and we've enabled learning
    break unless ($bad_starting_genes.include?(start_gene) and $learn_bad_starting_genes)
    # now if all genes are bad, well pick one randomly
    break if num_try > 10 * $top_limit
    # we have tried but found a bad gene
    num_try += 1
  end
  
  return start_gene
end

# --------------------------------------------------------------------------- #
# Cost function
# --------------------------------------------------------------------------- #
def cost_func cover, cov_pat_size
  w = 0
  miss = ($avail_pat.size - cov_pat_size).to_f/$avail_pat.size
  cover.each { |gene| w += $all_genes_weight[gene] }
  cost = $ALPHA*(miss) + $BETA*(w)
  
  return cost
end

# --------------------------------------------------------------------------- #
# Pick and Randomize Patients
# --------------------------------------------------------------------------- #

# randomly select the test set in the three part train/validation/test set
def pick_test_set
  $all_avail_genes = deep_copy_hash $all_genes
  $all_avail_pat   = deep_copy_hash $all_pat
  
  # randomly withhold 10% of the patients
  $HIDE_PERCENT = 0.1
  pick_patients
  
  # use the remaining 90% for range alpha runs
  $all_avail_genes = deep_copy_hash $genes
  $all_avail_pat   = deep_copy_hash $avail_pat
  $test_pat        = deep_copy_hash $hide_pat
  
  # next hide 20% for validation set
  $HIDE_PERCENT = 0.2
end

# randomly select % of the patients
def pick_patients
  $genes = Hash.new {|hash, key| hash[key] = Set.new() }
  $hide_pat = Hash.new {|hash, key| hash[key] = Set.new() }
  $avail_pat = Hash.new {|hash, key| hash[key] = Set.new() }

  pick_pat = Set.new($all_avail_pat.keys.sample(($HIDE_PERCENT*$all_avail_pat.size).to_i))
  $all_avail_pat.each_key do |pat|
    if pick_pat.include? pat
      $hide_pat[pat] = $all_avail_pat[pat]
    else
      $avail_pat[pat] = $all_avail_pat[pat]
      $avail_pat[pat].each { |gene| $genes[gene] << pat }
    end
  end
end

# --------------------------------------------------------------------------- #
# Execute nCOP heuristic functions
# --------------------------------------------------------------------------- #

# generic entry function for all algorithms
def exec_search alg
  cover = send("exec_search_#{alg}") 
  write_cover(cover, alg) if $extensive_write
  return cover
end

# driver for nCOP heruistic search
def exec_search_nCOP_heruistic
  # select a starting node
  s_gene = get_start_gene
  
  # initialize cover and other variables
  cover, cov_pat, $explored = [], Set.new, Set.new
  cover << s_gene and cov_pat.merge $genes[s_gene]
  start_cost = cost_func cover, cov_pat.size
  
  # run nCOP greedy heuristic
  nCOP_heuristic(cover, cov_pat, start_cost)
  
  # learn which genes are bad starting nodes
  if cover.size < 3
    $bad_starting_genes << s_gene
  end
  
  # return the cover found
  cover
end

# --------------------------------------------------------------------------- #
# nCOP Heuristic
# --------------------------------------------------------------------------- #

def nCOP_heuristic cover, cov_pat, old_cost
  # if we have covered all patients then we are done
  if cov_pat.size == $avail_pat.size then
    return
  end
  
  # Try expanding in one step
  min, max_g = $MIN_c, []

  cover.each do |v|
    $g.neighbors(v).each do |w|
      next if cover.include? w
        
      # compute the new cost and delta
      new_cost = cost_func Set.new(cover) << w, ((newly_cov_pat(w, cov_pat)).size + cov_pat.size)
      delta = new_cost - old_cost
      
      # if a better value, remember it
      if delta < min
        max_g = ([] << w)
        min = delta
      
      # if a tie, add it to the current best
      elsif delta == min
        max_g << w
      end
    end
  end
  
  # if we can expand then do so
  if min < $MIN_c then
    add_gene = max_g[rand(max_g.size)]
    cover <<  add_gene and cov_pat.merge((newly_cov_pat(add_gene, cov_pat)))

    # recursive call
    nCOP_heuristic(cover, cov_pat, cost_func(cover, cov_pat.size))
    return
  end
  
  # Try expanding in two steps
  max_g = []
  cover.each do |v|
    $g.neighbors(v).each do |w|
      next if cover.include? w
      
      $g.neighbors(w).each do |gene|
        next if cover.include? gene or !$genes.include? gene
   
        # compute the new cost and delta
        new_cost = cost_func(Array.new(cover) << w << gene, ((newly_cov_pat(gene, cov_pat).merge(newly_cov_pat(w, cov_pat))).size + cov_pat.size))
        delta = new_cost - old_cost
        pair = ([] << w << gene)
      
        # if a better value, remember it
        if delta < min
          max_g = ([] << pair)
          min = delta
          
        # if a tie, add it to the current best
        elsif delta == min
          max_g << pair
        end
      end
    end
  end
  
  # if we can expand do so
  if min < $MIN_c then
    add_gene = max_g[rand(max_g.size)]
    cover << add_gene[0] << add_gene[1] and cov_pat.merge(newly_cov_pat(add_gene[1], cov_pat).merge(newly_cov_pat(add_gene[0], cov_pat)))
    
    # recursive call
    nCOP_heuristic(cover, cov_pat, cost_func(cover, cov_pat.size))
    return
  end
  
  # we were not able to expand, terminate the recursion
  return
end

# --------------------------------------------------------------------------- #
# Execute Range Alpha
# --------------------------------------------------------------------------- #

def run_range_alpha
  $NUM_AVG_ITER = 100
  $alpha_range = [0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9]
  $run_alpha_range = true
  
  pick_test_set
  exec_range_alpha
end

# vary alpha and beta
def exec_range_alpha
  if $run_alpha_range
    puts "Iterating over range of alpha values. This may take a while..."
  else
    puts "Running nCOP with alpha = #{$alpha_range[0]}. Prioritizing genes..."
  end
  
  aggrStats = RunInfo.new
  
  # for each alpha and beta
  $alpha_range.each do |n|
    start_alp = Time.now
    
    # scale alpha and beta appropriately
    $ALPHA = n.round(4)
    $BETA = ((1 - $ALPHA).to_f/$normalization_const).round(4)
    
    # initialize an empty array for bad starting nodes
    $bad_starting_genes = []

    # run it NUM_ITER times
    exec_iter_fixed_alpha aggrStats
    
    puts "  alpha = #{n} done"
    #puts "avg iter for $ALPHA = #{n} completed after = #{proper_time_since(start_alp)}" if $run_alpha_range
  end
  
  # remember the stats
  $aggrStats = aggrStats
  
  puts "Done."
end

# for a given alpha runs NUM_ITER the algos and averages the values
def exec_iter_fixed_alpha aggrStats
  iterStats = RunInfo.new
  iterConsensus = Hash.new {|hash, key| hash[key] = 0 }
  
  # NUMBER_OF_ITER times do
  $NUM_AVG_ITER.times do |i| 
    # randomly withhold 20% of data
    pick_patients
    
    # find G'
    cover = exec_search "nCOP_heruistic"
      
    # store the stats
    iterStats.fillIn cover

    # add the cover to the consensus list
    add_to_consensus(iterConsensus, cover)
  end
  
  # compute the average and std for this set of iterations
  aggrStats.fillAggregateStats iterStats
  aggrStats.alpha << $ALPHA
  aggrStats.beta << $BETA

  # remember the frequency genes occured
  $selected_genes = iterConsensus
end

# --------------------------------------------------------------------------- #
# Select Alpha Procedure
# --------------------------------------------------------------------------- #

def pick_alpha aggrStats = $aggrStats
  puts ""
  puts "Selecting optimal alpha:"
  pick_alpha_rules aggrStats
  $alpha = sprintf('%.3f', $alpha_pick[$cancer])
  puts "$alpha = #{$alpha}"
  puts ""
end


def pick_alpha_rules aggrStats
  $alpha_pick  = {}
  alpha_chosen = aggrStats.alpha.map {|e| e.to_f.round(2)}
  hsh_hide     = aggrStats.hide_avg.map {|e| e.to_f.round(4)}
  hsh_aval     = aggrStats.aval_avg.map {|e| e.to_f.round(4)}
  
  # Rule Constants
  $len_inerval = 3
  $PLATO_WIDTH = 0.05
  $ALPHA_MAX_DIFF = 0.10
  
  # the alpha for which hide_avg is maximum
  alpha_max = alpha_chosen[hsh_hide.index(hsh_hide.max)]

  # find an alpha point
  (0..(alpha_chosen.size-$len_inerval)).each do |i|    
    plato = true
    
    # check if we are in a plato
    (1..$len_inerval).each do |j|
      if (hsh_hide[i].to_f - hsh_hide[i+j].to_f).abs > $PLATO_WIDTH
        plato = false
      end
    end
    
    if plato
      # if you are the first point in a plato with overfitting > large return you
      if (hsh_hide[i].to_f - hsh_aval[i].to_f).abs > $PLATO_WIDTH
        # unless you are more than large below max alpha, in which case keep going
        next if (hsh_hide[i].to_f - hsh_hide.max).abs > $ALPHA_MAX_DIFF
        
        $alpha_pick[$cancer] = alpha_chosen[i].to_s
        return
        
      # else if the next point is overfitting still return you
      elsif (hsh_hide[i+1].to_f - hsh_aval[i+1].to_f).abs > $PLATO_WIDTH
        # unless you are more than large below max alpha, in which case keep going
        next if (hsh_hide[i].to_f - hsh_hide.max).abs > $ALPHA_MAX_DIFF
        
        $alpha_pick[$cancer] = alpha_chosen[i].to_s
        return
        
      # else keep going
      else
        next
      end
    end
  end
  
  # well if none is good, select the highest point then, unless it is at the end
  if alpha_max < 0.84
    $alpha_pick[$cancer] = alpha_max.to_s
  else
    $alpha_pick[$cancer] = "0.8"
  end
end


# --------------------------------------------------------------------------- #
# Run nCOP given a selected alpha to prioritize genes
# --------------------------------------------------------------------------- #

def run_nCOP_given_alpha selected_alpha = $alpha
  # use the given alpha
  $alpha_range = ([] << selected_alpha.to_f.round(2))
  
  # set the rest of the params
  $NUM_AVG_ITER = 1000
  $run_alpha_range = false
  
  $all_avail_genes = deep_copy_hash $all_genes
  $all_avail_pat   = deep_copy_hash $all_pat
  
  # randomly withhold 15% of the patients
  $HIDE_PERCENT = 0.15
  
  exec_range_alpha
  
  # write the genes sorted by their frequency
  puts "Writing output to: #{$OUT_DIR}#{$cancer}_results.txt"
  genes = sort_hash($selected_genes)
  fo = File.open("#{$OUT_DIR}#{$cancer}_results.txt", 'w')
  genes.each { |e| fo.puts "#{e.ljust(10)} #{($selected_genes[e].to_f*100/$NUM_AVG_ITER).round(1)}%" }
  puts "Done."
end


# --------------------------------------------------------------------------- #
# RunInfo Class
# --------------------------------------------------------------------------- #
class RunInfo
  attr_accessor :aval_avg
  attr_accessor :aval_std
  attr_accessor :hide_avg
  attr_accessor :hide_std
  attr_accessor :test_avg
  attr_accessor :test_std
  attr_accessor :node_avg
  attr_accessor :node_std
  attr_accessor :alpha
  attr_accessor :beta
  
  def initialize
    @aval_avg, @aval_std, @hide_avg, @hide_std = [], [], [], []
    @test_avg, @test_std, @node_avg, @node_std = [], [], [], []
    @alpha, @beta = [],[]
  end
  
  def fillIn cover
    @aval_avg << get_frac_covered($avail_pat, cover)
    @hide_avg << get_frac_covered($hide_pat, cover)
    @test_avg << get_frac_covered($test_pat, cover)
    @node_avg << cover.size
  end
  
  def fillAggregateStats r
    @aval_avg << r.aval_avg.avg
    @aval_std << r.aval_avg.std
    @hide_avg << r.hide_avg.avg
    @hide_std << r.hide_avg.std
    @test_avg << r.test_avg.avg
    @test_std << r.test_avg.std
    @node_avg << r.node_avg.avg
    @node_std << r.node_avg.std
  end
  
  def write_runInfo alg = "nCOP"
    fo = File.open("#{$OUT_DIR}#{$cancer}_AVG.txt", 'w')
    fo.puts \
    "#{alg}_aval_avg = #{@aval_avg.inspect}",
    "#{alg}_aval_std = #{@aval_std.inspect}",
    "#{alg}_hide_avg = #{@hide_avg.inspect}",
    "#{alg}_hide_std = #{@hide_std.inspect}",
    "#{alg}_test_avg = #{@test_avg.inspect}",
    "#{alg}_test_std = #{@test_std.inspect}",
    "#{alg}_node_avg = #{@node_avg.inspect}",
    "#{alg}_node_std = #{@node_std.inspect}",
    "alphas = #{@alpha.inspect}",
    "betas  = #{@beta.inspect}",
    ""
    fo.close
  end
end

# add the selected genes to the consensus
def add_to_consensus iterConsensus, cover
  cover.each { |e| iterConsensus[e] += 1 }
end


# --------------------------------------------------------------------------- #
# Small helper functions
# --------------------------------------------------------------------------- #

# return the fraction of covered patients
def get_frac_covered patients, cover
  return 0 if patients.size == 0
  (patients.keys.count { |pat| !(patients[pat]&cover).empty? }.to_f/patients.size).round($d_prec)
end

# return the number of covered patients
def get_num_covered patients, cover
  patients.keys.count { |pat| !(patients[pat]&cover).empty? }
end

# return the number of newly covered patients by the specified gene
def newly_cov_pat gene, cov_pat
  new_cov = Set.new
  return new_cov unless $genes.has_key? gene
  $genes[gene].each { |pat| new_cov << pat unless cov_pat.include? pat }
  new_cov 
end

# make a deep copy of a hash
def deep_copy_hash g
  h = {}
  g.each_pair do |k,v|
    v = Marshal.load(Marshal.dump(v))
    h[k] = v
  end
  h
end

# sort a hash
def sort_hash hsh
  arr = []
  hsh.values.uniq.sort.reverse.each do |v|
    hsh.map{ |k,v1| (v1==v) ? k : nil }.compact.each { |e| arr << e }
  end
  arr
end

# return nicely formatted time since the start time
def proper_time_since start_time
  total_time = Time.now.to_f - start_time.to_f
  return "#{total_time.round(2)} sec" if total_time < 60
  return "#{(total_time.to_f/60).round(2)} min" if total_time < 3600
  return "#{(total_time.to_f/3600).round(2)} h"
end


# --------------------------------------------------------------------------- #
# Array Class
# --------------------------------------------------------------------------- #
class Array
  def sum
    self.inject(:+)
  end
  
  def avg
    (self.inject(:+).to_f/self.size).round($d_prec)
  end
  
  def std
    mu = self.avg
    n = self.size - (self.size > 25 ? 1 : 0)
    (Math.sqrt(self.map { |r| (r - mu) ** 2}.inject(:+).to_f/n)).round($d_prec)
  end
  
  def median
      sorted = self.sort
      len = sorted.length
      return (sorted[(len - 1) / 2] + sorted[len / 2]) / 2.0
  end
  
  def avg_column i
    self.map { |r| r[i]}.inject(:+).to_f/self.size
  end
  
  def std_column i
    mu = self.avg_column i
    n = self.size - (self.size > 25 ? 1 : 0)
    Math.sqrt(self.map { |r| (r[i]- mu) ** 2}.inject(:+).to_f/n)
  end

  def ratio_col i, j
    self.map { |r| r[i].to_f/r[j]}
  end
end

# --------------------------------------------------------------------------- #
# Simple Graph Class
# --------------------------------------------------------------------------- #
class RubGraph
  attr_reader :nodes
  
  def initialize(id = 'id_default')
    @id = id
    @nodes = {}
  end

  def add_node n
    @nodes[n] = {} unless @nodes.include? n
  end
  
  def add_edge e1, e2
    add_node e1; add_node e2
    @nodes[e1][e2] = 1 and @nodes[e2][e1] = 1
  end
  
  def remove_node n
    raise "node #{n} not in the graph" unless @nodes.has_key?(n)
    @nodes.delete n
    @nodes.each { |_, v| v.delete n }
  end
  
  def remove_edge n1, n2
    @nodes[n1].delete n2
    @nodes[n2].delete n1
  end
  
  def size
    @nodes.size
  end
  
  def degree v
    @nodes[v].size
  end
  
  def neighbors v
    @nodes[v].keys
  end
  
  def all_nodes_except v
    @nodes.keys.clone - ([]<<v)
  end
  
  def num_nodes
    @nodes.keys.size
  end
  
  def num_edges
    sum = 0
    @nodes.keys.each do |v|
      sum += @nodes[v].keys.size
    end
    sum/2
  end
end

# --------------------------------------------------------------------------- #
# The Pickup Library
# --------------------------------------------------------------------------- #

class Pickup
  attr_reader :list, :uniq
  attr_writer :pick_func, :key_func, :weight_func

  def initialize(list, opts={}, &block)
    @list = list
    @uniq = opts[:uniq] || false
    @pick_func = block if block_given?
    @key_func = opts[:key_func]
    @weight_func = opts[:weight_func]
  end

  def pick(count=1, opts={}, &block)
    func = block || pick_func
    key_func = opts[:key_func] || @key_func
    weight_func = opts[:weight_func] || @weight_func
    mlist = MappedList.new(list, func, uniq: uniq, key_func: key_func, weight_func: weight_func)
    result = mlist.random(count)
    count == 1 ? result.first : result
  end

  def pick_func
    @pick_func ||= begin
      Proc.new do |val|
        val
      end
    end
  end

  class CircleIterator
    attr_reader :func, :obj, :max, :key_func, :weight_func

    def initialize(obj, func, max, opts={})
      @obj = obj.dup
      @func = func
      @max = max
      @key_func = opts[:key_func] || key_func
      @weight_func = opts[:weight_func] || weight_func
    end

    def key_func
      @key_func ||= begin
        Proc.new do |item|
          item[0]
        end
      end
    end

    def weight_func
      @weight_func ||= begin
        Proc.new do |item|
          item[1]
        end
      end
    end

    def each
      until obj.empty?
        start = 0
        obj.each do |item|
          key = key_func.call(item)
          weight = weight_func.call(item)

          val = func.call(weight)
          start += val
          if yield([key, start, max])
            obj.delete key
            @max -= val
          end
        end
      end
    end
  end

  class MappedList
    attr_reader :list, :func, :uniq, :key_func, :weight_func

    def initialize(list, func, opts=nil)
      if Hash === opts
        @key_func = opts[:key_func]
        @weight_func = opts[:weight_func] || weight_func
        @uniq = opts[:uniq] || false
      else
        if !!opts == opts
          # If opts is explicitly provided as a boolean, show the deprecated warning.
          warn "[DEPRECATED] Passing uniq as a boolean to MappedList's initialize method is deprecated. Please use the opts hash instead."
        end

        @uniq = opts || false
      end

      @func = func
      @list = list
      @current_state = 0
    end

    def weight_func
      @weight_func ||= begin
        Proc.new do |item|
          item[1]
        end
      end
    end

    def each(&blk)
      CircleIterator.new(@list, func, max, key_func: @key_func, weight_func: weight_func).each do |item|
        if uniq
          true if yield item
        else
          nil while yield(item)
        end
      end
    end

    def random(count)
      raise "List is shorter then count of items you want to get" if uniq && list.size < count
      nums = count.times.map{ rand(max) }.sort
      return [] if max == 0
      get_random_items(nums)
    end

    def get_random_items(nums)
      current_num = nums.shift
      items = []
      each do |item, counter, mx|
        break unless current_num
        if counter%(mx+1) > current_num%mx
          items << item
          current_num = nums.shift
          true
        end
      end
      items
    end

    def max
      @max ||= begin
        max = 0
        list.each{ |item| max += func.call(weight_func.call(item)) }
        max
      end
    end
  end
end

# --------------------------------------------------------------------------- #
# Required libraries and basic paths and constants 
# --------------------------------------------------------------------------- #
require 'set'
require 'fileutils'

# Set the base directory to that containing this file
$BASE_DIR = (File.dirname $0) + "/"
Dir.chdir $BASE_DIR

# Set the output directory to BASE_DIR/Outputs
$OUT_DIR = "Outputs/"

$MIN_c, $d_prec = -0.000_001, 4
$top_limit = 5
$learn_bad_starting_genes = true
$picked = Hash.new {|hash, key| hash[key] = 0 }
$bad_starting_genes = []
$HIDE_PERCENT = 0.2


if __FILE__ == $PROGRAM_NAME
  read_user_input
  run_according_to_user_choices
end
