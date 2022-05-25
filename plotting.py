import plotly.graph_objects as go

COLORS=['rgb(31, 119, 180)', 'rgb(255, 127, 14)',
        'rgb(44, 160, 44)', 'rgb(214, 39, 40)',
        'rgb(148, 103, 189)', 'rgb(140, 86, 75)',
        'rgb(227, 119, 194)', 'rgb(127, 127, 127)',
        'rgb(188, 189, 34)', 'rgb(23, 190, 207)']

def cycle(my_list):
    start_at = 0
    while True:
        yield my_list[start_at]
        start_at = (start_at + 1) % len(my_list)

class Plot:
    def extracted_chroms(self, df_chrom, chroms=[], df_auc = None):
        colors = cycle(COLORS)
        traces_line = []
        traces_bar = []
        for chrom in chroms:
            color = next(colors)
            traces_line.append(go.Scatter(x=df_chrom["time"], y=df_chrom[chrom], name=chrom, mode="lines", line_color=color))
            if len(df_auc) == 1 and chrom in df_auc.columns:
                traces_bar.append(go.Bar(x=[chrom], y=[df_auc[chrom][0]], name=chrom, marker_color=color))
        fig_chrom = go.Figure()
        fig_auc = go.Figure()
        for trace in traces_line:
            fig_chrom.add_trace(trace)
        for trace in traces_bar:
            fig_auc.add_trace(trace)
        fig_chrom.update_layout(xaxis=dict(title="time"), yaxis=dict(title="intensity (cps)"))  
        fig_auc.update_layout(xaxis=dict(title=""), yaxis=dict(title="area under curve (counts)"))
        fig_auc.update_traces(width=0.15)
        return fig_chrom, fig_auc      

    def FFMID_chroms_from_df(self, df, compounds = []):
        colors = cycle(COLORS)
        if compounds:
            trace_names = compounds
        else:
            trace_names = set([c[:-5] for c in df.columns if c.endswith("RT")])
        traces = []
        for name in trace_names:
            color = next(colors)
            i = 1
            while name+"_"+str(i)+"_RT" in df.columns:
                label = name+"_"+str(i)
                traces.append(go.Scatter(x=df[label+"_RT"], y=df[label+"_int"], name=label, mode="lines", line_color=color)) 
                i += 1
        fig = go.Figure()
        for trace in traces:
            fig.add_trace(trace)
        fig.update_layout(xaxis=dict(title="time"), yaxis=dict(title="intensity (cps)"))
        return fig