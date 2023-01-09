import pandas as pd
import numpy as np

class Statistics:
    def normalize_max(self, df):
        """Normalize values between 0 and 1 for the complete DataFrame.
        Parameters
        ----------
        df : pandas.DataFrame
            Input DataFrame with non normalized values.
        Returns
        -------
        pandas.DataFrame
            Normalized DataFrame (based on the maximum value in a DataFrame).
        """
        column_maxes = df.max()

        df_max = column_maxes.max()

        return df / df_max


    def maximum_absolute_scaling_per_column(self, df):
        """Normalize values between 0 and 1 for each column seperately.
        Parameters
        ----------
        df : pandas.DataFrame
            Input DataFrame with non normalized values.
        Returns
        -------
        pandas.DataFrame
            Normalized DataFrame (each column is normalized on it's maximum value)
        """
        # copy the dataframe
        df_scaled = df.copy()
        # apply maximum absolute scaling
        for column in df_scaled.columns:
            column_values = pd.to_numeric(df_scaled[column])
            df_scaled[column] = column_values / column_values.abs().max()
        return df_scaled


    def get_mean_std_change_df(self, df, sample_pairs):
        """Calculate DataFrames for mean, standard deviation and fold change.
        Given a DataFrame and a list of sample pairs (that will be compared in fold change)
        DataFrames for mean values between replicates, their standard deviation and fold changes
        are calculated.
        Parameters
        ----------
        df : pandas.DataFrame
            Input DataFrame with ID names as index and sample names (without the replicate specification)
            as columns. Only columns specified in sample_pairs will be considered for calculation!
        sample_pairs : list of tuples
            Sample names to be compared with each other for fold change calculations (eg. [('control', 'treatment')])
            will calculate the log2 fold change of normalized (on maximum) values for the means of treatment / control.
        Returns
        -------
        df_mean : pandas.DataFrame
            DataFrame with mean values from the replicates per sample.
        df_std : pandas.DataFrame
            DataFrame with standard deviation values from the replicates per sample.
        df_change : pandas.DataFrame
            DataFrame with log2 fold change values for the specified sample pairs.
        """

        df_mean = pd.DataFrame(index=df.index)
        df_std = pd.DataFrame(index=df.index)
        df_change = pd.DataFrame(index=df.index)
        
        # replicate numbers should be at the end of the file name separated with # (e.g. samplename#1)
        for name in [sample for sample_pair in sample_pairs for sample in sample_pair]:
            replicates = [c for c in df.columns if c.startswith(name + "#") or c == name]
            print(replicates)
            df_mean[name] = df[replicates].mean(axis=1)
            df_std[name] = df[replicates].std(axis=1)

        for pair in sample_pairs:
            df_pair = self.normalize_max(df_mean[[pair[0], pair[1]]] + 1)
            df_change[pair[1] + "/" + pair[0]] = np.log2(
                df_pair[pair[1]] / df_pair[pair[0]]
            )
        df_mean = df_mean.applymap(lambda x: int(round(x, 0)) if isinstance(x, (int, float, )) else x)
        df_std.fillna(0, inplace=True)
        df_std = df_std.applymap(lambda x: int(round(x, 0)) if isinstance(x, (int, float)) else x)
        return df_mean, df_std, df_change
