import plotly.graph_objects as go
import plotly.express as px
import pandas as pd

COLORS=['#636efa', '#ef553b',
        '#00cc96', '#ab63fa',
        '#ffa15a', '#19d3f3',
        '#ff6692', '#b6e880',
        '#fe9aff', 'rgb(23, 190, 207)']

def cycle(my_list):
    start_at = 0
    while True:
        yield my_list[start_at]
        start_at = (start_at + 1) % len(my_list)

class Plot:
    def extracted_chroms(self, df_chrom, chroms=[], df_auc = None, title="", time_unit="seconds"):
        colors = cycle(COLORS)
        traces_line = []
        traces_bar = []
        for chrom in chroms:
            if chrom == "BPC":
                color = "#CCCCCC"
            elif chrom == "AUC baseline":
                color = "#555555"
            else:
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
        fig_chrom.update_layout(title=title, xaxis=dict(title=f"time ({time_unit})"), yaxis=dict(title="intensity (counts per second)"))  
        fig_auc.update_layout(title=title, xaxis=dict(title=""), yaxis=dict(title="area under curve (counts)"))
        fig_auc.update_traces(width=0.3)
        return fig_chrom, fig_auc        

    def FFMID(self, df_chrom, compounds = [], df_auc = None, df_auc_combined = None, title="", time_unit="seconds"):
        colors = cycle(COLORS)
        if compounds:
            trace_names = compounds
        else:
            trace_names = set([c[:-5] for c in df_chrom.columns if c.endswith("RT")])
        traces_line = []
        traces_bar = []
        traces_bar_combined = []
        for name in trace_names:
            color = next(colors)
            i = 1
            while name+"_"+str(i)+"_RT" in df_chrom.columns:
                label = name+"_"+str(i)
                traces_line.append(go.Scatter(x=df_chrom[label+"_RT"], y=df_chrom[label+"_int"], name=label, mode="lines", line_color=color)) 
                if len(df_auc) == 1 and name in df_auc.columns and i == 1:
                    traces_bar.append(go.Bar(x=[name], y=[df_auc[name][0]], name=name, marker_color=color))
                i += 1
        if len(df_auc) == 1:
            colors = cycle(COLORS)
            for name in df_auc_combined.columns:
                color = next(colors)
                traces_bar_combined.append(go.Bar(x=[name], y=[df_auc_combined[name][0]], name=name, marker_color=color))
        fig_chrom = go.Figure()
        fig_auc = go.Figure()
        fig_auc_combined = go.Figure()
        for trace in traces_line:
            fig_chrom.add_trace(trace)
        for trace in traces_bar:
            fig_auc.add_trace(trace)
        for trace in traces_bar_combined:
            fig_auc_combined.add_trace(trace)
        fig_chrom.update_layout(title=title, xaxis=dict(title=f"time ({time_unit})"), yaxis=dict(title="intensity (cps)"))
        for fig in (fig_auc, fig_auc_combined):
            fig.update_layout(title=title, xaxis=dict(title=""), yaxis=dict(title="area under curve (counts)"))
            fig.update_traces(width=0.3)
        return fig_chrom, fig_auc, fig_auc_combined

    def FeatureMatrix(self, df, df_std=pd.DataFrame(), samples=[], features=[], title="", y_title = "area under curve (counts)"):
        samples = df.columns.tolist()
        features = df.index.tolist()
        fig = go.Figure()
        for feature in features:
            if not df_std.empty:
                fig.add_trace(go.Bar(x=samples, y=df.loc[feature], name=feature, error_y=dict(type="data", array=df_std.loc[feature], visible=True)))
            else:
                fig.add_trace(go.Bar(x=samples, y=df.loc[feature], name=feature))
        fig.update_layout(title=title, yaxis=dict(title=y_title))
        return fig
    
    def FeatureMatrixHeatMap(self, df, title=""):
        fig = go.Figure(data=go.Heatmap({'z': df.values.tolist(),
                                        'x': df.columns.tolist(),
                                        'y': df.index.tolist()},
                                        colorscale= [[0, 'rgba(69, 117, 180, 1.0)'],   
                                               [ 0-df.min().min()/(df.max().max()-df.min().min()), 'rgba(255, 255, 255, 1.0)'],  
                                        [1, 'rgba(215,48,39, 1.0)']]))
        # fig.layout.width = 200+10*len(df.columns)
        # fig.layout.height = 30*len(df.index)
        fig.update_layout(title=title)
        return fig
