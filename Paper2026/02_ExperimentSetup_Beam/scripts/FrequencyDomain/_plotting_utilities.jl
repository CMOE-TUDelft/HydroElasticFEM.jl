using Plots

# ============================================================================
# Base minimal plotting utilities for Frequency Domain simulations
# ============================================================================

function plot_power_coefficient(cache_for_plots)
  
  @info """
  Plotting preliminary power coefficient\n
  Fields in cache_for_plots:
  $( fieldnames(typeof(cache_for_plots)) ) \n
  """

  @unpack ω, ω_nat, prbPow, filename = cache_for_plots
  
  @show filename

  mask = map(x -> x in ω_nat, ω)
  idx = findall(mask)

  # Power Balance Plots
  plt1 = plot(ω, prbPow[:,3]./prbPow[:,2], linewidth=3, 
    xlabel = "ω (rad/s)",
    ylabel = "K_R",
    title = "Reflection coefficient",
    ylims = (0,1.0))
  
  plot!(plt1, ω[idx], prbPow[idx,3]./prbPow[idx,2],
    seriestype = :scatter, markersize=3,
    label = "ω_nat")

  plt2 = plot(ω, prbPow[:,4]./prbPow[:,2], linewidth=3, 
    xlabel = "ω (rad/s)",
    ylabel = "K_T",
    title = "Transmission coefficient",
    ylims = (0,1.0))

  plot!(plt2, ω[idx], prbPow[idx,4]./prbPow[idx,2],
    seriestype = :scatter, markersize=3,
    label = "ω_nat")
  

  plt3 = plot(ω, prbPow[:,5]./prbPow[:,2], linewidth=3, 
    xlabel = "ω (rad/s)",
    ylabel = "K_A",
    title = "Absorption coefficient",
    ylims = (0,1.0))
  
  plot!(plt3, ω[idx], prbPow[idx,5]./prbPow[idx,2],
    seriestype = :scatter, markersize=3,
    label = "ω_nat")

  plt4 = plot(ω, 100 * prbPow[:,6]./prbPow[:,2], linewidth=3, 
    xlabel = "ω (rad/s)",
    ylabel = "Error %",
    title = "Power Relative Error")
  
  plot!(plt4, ω[idx], prbPow[idx,6]./prbPow[idx,2],
    seriestype = :scatter, markersize=3,
    label = "ω_nat")

  
  pltAll = plot(plt1, plt2, plt3, plt4, layout=4, dpi=330,
    plot_title = "Power Balance")
  savefig(pltAll,filename*"_plots/mem_powerBalance"*".png")

  println()
  return 0
end


function plot_resonator_RAO(cache_for_plots)
  
  @info """
  Plotting resonator RAO\n
  Fields in cache_for_plots:
  $( fieldnames(typeof(cache_for_plots)) ) \n
  """

  @unpack ω, ω_nat, prbResnRAO, filename = cache_for_plots
  
  @show filename

  mask = map(x -> x in ω_nat, ω)
  idx = findall(mask)

  # Resonator RAO Plots
  plt3 = plot(ω, abs.(prbResnRAO[:,2]), linewidth=3,
    xlabel = "ω (rad/s)",
    ylabel = "q_resn (m)",
    title = "Resonator Displacement Amplitude")
  
  plot!(plt3, ω[idx], abs.(prbResnRAO[idx,2]),
    seriestype = :scatter, markersize=3,
    label = "ω_nat")
  
  plt4 = plot(ω, abs.(prbResnRAO[:,3]), linewidth=3,
    xlabel = "ω (rad/s)",
    ylabel = "η_resn (m)",
    title = "Membrane Elevation at Resonator Location")
  
  plot!(plt4, ω[idx], abs.(prbResnRAO[idx,3]),
    seriestype = :scatter, markersize=3,
    label = "ω_nat")
  
  pltAll = plot(plt3, plt4, layout=(2,1), dpi=330,
    plot_title = "Resonator RAO")
  savefig(pltAll,filename*"_plots/mem_resonatorRAO"*".png")

  println()
  return 0
end


function plot_probes_along_free_surface(cache_for_plots)

  @info """
  Plotting probes along the free surface\n
  Fields in cache_for_plots:
  $( fieldnames(typeof(cache_for_plots)) ) \n
  """

  @unpack filename, ω, prbDa, prbDa_x, prbxy, prbx = cache_for_plots

  
  @show filename

  for lprb in 1:length(prbxy)    

    plt1 = plot(ω, abs.(prbDa[:,lprb]), linewidth=3, 
      xlabel = "ω (rad/s)",
      ylabel = "A (m)",
      title = "Amplitude")  

    plt2 = plot(ω, abs.(prbDa_x[:,lprb]), linewidth=3, 
      xlabel = "ω (rad/s)",
      ylabel = "dA/dx",
      title = "Slope Magnitude")
    
    plt3 = plot(ω, angle.(prbDa[:,lprb]), linewidth=3, 
      xlabel = "ω (rad/s)",
      ylabel = "α (rad)",
      title = "Phase")  

    plt4 = plot(ω, angle.(prbDa_x[:,lprb]), linewidth=3, 
      xlabel = "ω (rad/s)",
      ylabel = "α (rad)",
      title = "Slope Phase")
    
    xloc = prbx[lprb]
    pltAll = plot(plt1, plt2, plt3, plt4, layout=4, dpi=330,
      plot_title = "x = $xloc")

    savefig(pltAll,filename*"_plots/mem_dxPrb_$lprb"*".png")
  end

  println()
  return 0
end  