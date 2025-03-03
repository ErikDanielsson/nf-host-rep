params.analysis_dir = "analysis"

process generate_trace_plots {
    publishDir "analysis"

    input:
        path revbayes_out 
        path treeppl_out


    output:
        path "trace_plot.png"
    
    script:
    """
    trace_plots.py ${revbayes_out.join(' ')} ${treeppl_out.join(' ')}
    """
}