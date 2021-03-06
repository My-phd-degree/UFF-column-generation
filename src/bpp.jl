using CPLEX

mutable struct DataBPP
  items::Array{Float64}
  C::Float64
end

function solveBPP(dataBPP::DataBPP)
  M = Model(solver = CplexSolver(
#                                 CPX_PARAM_MIPDISPLAY=0,
#                                 CPX_PARAM_SCRIND=0
                                ))
  n = length(dataBPP.items)
  Bins = 1:n
  @variable(M, 0 <= y[j in Bins] <= 1, Int)
  @variable(M, 0 <= x[(i, j) in [(i, j) for i in 1:n for j in Bins]] <= 1, Int)
  @objective(M, Min, sum(y[j] for j in Bins))
  @constraint(M, every_item_must_be_packed_only_once[i in 1:n], sum(x[(i, j)] for j in Bins) == 1)
  @constraint(M, bin_limit[j in Bins], sum(x[(i, j)] * dataBPP.items[i] for i in 1:n) <= dataBPP.C)
  @constraint(M, open_bin[i in 1:n, j in Bins], x[(i, j)] <= y[j])
  solve(M)
  return getobjectivevalue(M)
end
