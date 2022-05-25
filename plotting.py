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