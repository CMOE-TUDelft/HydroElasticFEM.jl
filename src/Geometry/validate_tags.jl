using GridapGmsh

# ─────────────────────────────────────────────────────────────────────────────
# Gmsh tag validator
# ─────────────────────────────────────────────────────────────────────────────

"""
    validate_gmsh_tags(
        msh_file::String;
        required::Vector{String} = STANDARD_TAGS,
        dim::Int = 2,
    ) -> Bool

Load `msh_file` with GridapGmsh and check that every tag in `required`
appears as a physical-group name.

Prints a formatted report to stdout showing which tags are present and
which are missing, then returns `true` if all required tags are present,
`false` otherwise.

# Arguments
- `msh_file`  — path to the `.msh` file
- `required`  — list of tag names to check; defaults to [`STANDARD_TAGS`](@ref)
- `dim`       — spatial dimension (used for display only)

# Example

```julia
ok = validate_gmsh_tags("tank.msh")
ok = validate_gmsh_tags("tank.msh"; required=["fluid","free_surface"])
```
"""
function validate_gmsh_tags(
  msh_file::String;
  required::Vector{String} = STANDARD_TAGS,
  dim::Int = 2,
)
  if !isfile(msh_file)
    println("ERROR: File not found: \"$msh_file\"")
    return false
  end

  model  = GridapGmsh.GmshDiscreteModel(msh_file)
  labels = get_face_labeling(model)
  present_names = Set(labels.tag_to_name)

  # Filter out internal Gridap labels
  detected = filter(
    n -> !startswith(n, "tag_") && n != "interior" && n != "boundary",
    collect(present_names),
  )

  println("─────────────────────────────────────────────────────────")
  println("validate_gmsh_tags: \"$msh_file\"  (dim=$dim)")
  println("─────────────────────────────────────────────────────────")
  println("  Detected physical groups: ",
    isempty(detected) ? "(none)" : join(sort(detected), ", "))
  println()

  all_ok  = true
  for tag in required
    present = tag in present_names
    status  = present ? "✓  PRESENT" : "✗  MISSING"
    println("  $status   \"$tag\"")
    all_ok = all_ok && present
  end

  extra = sort(setdiff(detected, required))
  if !isempty(extra)
    println()
    println("  Extra (non-required) tags:")
    for tag in extra
      println("      \"$tag\"")
    end
  end

  println()
  if all_ok
    println("  Result: ALL required tags present ✓")
  else
    missing_list = filter(t -> !(t in present_names), required)
    println("  Result: MISSING tags — ", join(missing_list, ", "))
  end
  println("─────────────────────────────────────────────────────────")
  return all_ok
end
