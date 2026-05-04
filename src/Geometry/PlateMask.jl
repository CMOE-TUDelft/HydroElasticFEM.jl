"""
    get_plate_triangulation(Γ, xb₀, xb₁, yb₀, yb₁)

Split a top-surface triangulation `Γ` into plate and free-surface
sub-triangulations by coordinate mask.

A cell belongs to the plate only if all its nodes satisfy:
`xb₀ ≤ x ≤ xb₁` and `yb₀ ≤ y ≤ yb₁`.

Returns `(Γb, Γf, Λb)` where `Λb = Skeleton(Γb)`.
"""
function get_plate_triangulation(Γ, xb₀, xb₁, yb₀, yb₁)
  function is_plate(cell_nodes)
    minimum([(xb₀ <= n[1] <= xb₁) && (yb₀ <= n[2] <= yb₁)
             for n in cell_nodes])
  end

  xΓ = get_cell_coordinates(Γ)
  mask = lazy_map(is_plate, xΓ)
  Γb_idx = findall(mask)
  Γf_idx = findall(!, mask)

  Γb = Triangulation(Γ, Γb_idx)
  Γf = Triangulation(Γ, Γf_idx)
  Λb = Skeleton(Γb)

  return Γb, Γf, Λb
end
