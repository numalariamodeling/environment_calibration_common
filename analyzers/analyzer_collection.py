##### Import required packages #####
# standard packages
import os
import pandas as pd
from idmtools.entities import IAnalyzer
from idmtools.entities.simulation import Simulation


class MonthlyIncidenceAnalyzer(IAnalyzer):
    """
    Take monthly MalariaSummaryReport and pull out the PfPR, Cases, Severe Cases
    and Population for each agebins.
    """

    def __init__(
        self,
        sweep_variables=None,
        working_dir="./",
        start_year=2000,
        end_year=2001,
        filter_exists=False,
    ):
        super(MonthlyIncidenceAnalyzer, self).__init__(
            working_dir=working_dir,
            filenames=[f"output/MalariaSummaryReport_Monthly_incidence_{x}.json" for x in range(start_year, end_year)],
        )
        self.sweep_variables = sweep_variables or ["Run_Number"]
        self.start_year = start_year
        self.end_year = end_year
        self.filter_exists = filter_exists

    def filter(self, simulation: Simulation):
        if self.filter_exists:
            file = os.path.join(simulation.get_path(), self.filenames[0])
            return os.path.exists(file)
        else:
            return True

    def map(self, data, simulation: Simulation):
        adf = []
        fname = self.filenames[0]

        for year, fname in zip(range(self.start_year, self.end_year), self.filenames):
            age_bins = data[fname]["Metadata"]["Age Bins"]
            df = pd.DataFrame.from_dict(
                data[fname]["DataByTimeAndAgeBins"], orient="columns"
            )[:-1]
            df["month"] = [[x] * len(age_bins) for x in range(1, len(df) + 1)]
            df["agebin"] = [age_bins] * len(df)
            df.rename(
                columns={
                    "Annual Clinical Incidence by Age Bin": "Cases",
                    "Average Population by Age Bin": "Pop"
                },
                inplace=True,
            )
            df = df[["agebin", "month", "Cases", "Pop"]]
            df = df.explode(list(df.columns))
            df["Year"] = year
            adf = adf + [df]

        adf = pd.concat(adf)

        for sweep_var in self.sweep_variables:
            if sweep_var in simulation.tags.keys():
                adf[sweep_var] = simulation.tags[sweep_var]

        return adf

    def reduce(self, all_data):
        selected = [data for sim, data in all_data.items()]
        if len(selected) == 0:
            print("\nWarning: No data have been returned... Exiting...")
            return

        adf = pd.concat(selected).reset_index(drop=True)
        adf.to_csv(
            (os.path.join(self.working_dir, "ClinicalIncidence_monthly.csv")),
            index=False,
        )

class MonthlyPfPRAnalyzer(IAnalyzer):
    """
    Take monthly MalariaSummaryReport and pull out the PfPR, Cases, Severe Cases
    and Population for each agebins.
    """

    def __init__(
        self,
        sweep_variables=None,
        working_dir="./",
        start_year=2000,
        end_year=2001,
        grp = 'U5',
        filter_exists=False,
    ):
        super(MonthlyPfPRAnalyzer, self).__init__(
            working_dir=working_dir,
            filenames=[
                f"output/MalariaSummaryReport_Monthly_prevalence_{grp}_{x}.json"
                for x in range(start_year, end_year)
            ],
        )
        self.sweep_variables = sweep_variables or ["Run_Number"]
        self.start_year = start_year
        self.end_year = end_year
        self.filter_exists = filter_exists
        self.grp = grp

    def filter(self, simulation: Simulation):
        if self.filter_exists:
            file = os.path.join(simulation.get_path(), self.filenames[0])
            return os.path.exists(file)
        else:
            return True

    def map(self, data, simulation: Simulation):
        adf = []
        fname = self.filenames[0]

        for year, fname in zip(range(self.start_year, self.end_year), self.filenames):
            age_bins = data[fname]["Metadata"]["Age Bins"]
            df = pd.DataFrame.from_dict(
                data[fname]["DataByTimeAndAgeBins"], orient="columns"
            )[:-1]
            df["month"] = [[x] * len(age_bins) for x in range(1, len(df) + 1)]
            df["agebin"] = [age_bins] * len(df)
            df.rename(
                columns={
                    "PfPR by Age Bin": "PfPR",
                    "Average Population by Age Bin": "Pop"
                },
                inplace=True,
            )
            df = df[["agebin", "month", "PfPR", "Pop"]]
            df = df.explode(list(df.columns))
            df["Year"] = year
            adf = adf + [df]

        adf = pd.concat(adf)

        for sweep_var in self.sweep_variables:
            if sweep_var in simulation.tags.keys():
                adf[sweep_var] = simulation.tags[sweep_var]

        return adf

    def reduce(self, all_data):
        selected = [data for sim, data in all_data.items()]
        if len(selected) == 0:
            print("\nWarning: No data have been returned... Exiting...")
            return

        adf = pd.concat(selected).reset_index(drop=True)
        adf.to_csv(
            (os.path.join(self.working_dir, f"PfPR_monthly_{self.grp}.csv")),
            index=False,
        )



class AnnualPfPRAnalyzer(IAnalyzer):
    """
    Take monthly MalariaSummaryReport and pull out the PfPR, Cases, Severe Cases
    and Population for each agebins.
    """

    def __init__(
        self,
        sweep_variables=None,
        working_dir="./",
        start_year=2000,
        end_year=2001,
        filter_exists=False,
    ):
        super(AnnualPfPRAnalyzer, self).__init__(
            working_dir=working_dir,
            filenames=[f"output/MalariaSummaryReport_Yearly_prevalence_{start_year}_to_{end_year}.json"]
        )

        self.sweep_variables = sweep_variables or ["Run_Number"]
        self.start_year = start_year
        self.end_year = end_year
        self.filter_exists = filter_exists

    def filter(self, simulation: Simulation):
        if self.filter_exists:
            file = os.path.join(simulation.get_path(), self.filenames[0])
            return os.path.exists(file)
        else:
            return True

    def map(self, data, simulation: Simulation):
        adf = []
        
        fname = self.filenames[0]

        for fname in self.filenames:
            
            age_bins = data[fname]["Metadata"]["Age Bins"]
            years = data[fname]["DataByTime"]["Time Of Report"]
            years = [y / 365 - 1 for y in years]
            
            df = pd.DataFrame.from_dict(data[fname]["DataByTimeAndAgeBins"],orient="columns")[:-1]
            #df["month"] = [[x] * len(age_bins) for x in range(1, len(df) + 1)]
            df["agebin"] = [age_bins] * len(df)
            df["Year"] = [[years[x]] * len(age_bins) for x in range(len(df))]
            df.rename(
               columns={
                   "PfPR by Age Bin": "PfPR",
                   "Average Population by Age Bin": "Pop"
               },
                inplace=True,
            )
            df = df[["agebin", "Year","PfPR", "Pop"]]
            df = df.explode(list(df.columns))
           
            adf = adf + [df]

            adf = pd.concat(adf)

        for sweep_var in self.sweep_variables:
            if sweep_var in simulation.tags.keys():
                adf[sweep_var] = simulation.tags[sweep_var]

        return adf

    def reduce(self, all_data):
        selected = [data for sim, data in all_data.items()]
        if len(selected) == 0:
            print("\nWarning: No data have been returned... Exiting...")
            return

        adf = pd.concat(selected).reset_index(drop=True)
        adf.to_csv(
            (os.path.join(self.working_dir, "PfPR_annual.csv")),
            index=False,
        )


class AnnualIncidenceAnalyzer(IAnalyzer):
    """
    Take monthly MalariaSummaryReport and pull out the PfPR, Cases, Severe Cases
    and Population for each agebins.
    """

    def __init__(
        self,
        sweep_variables=None,
        working_dir="./",
        start_year=2000,
        end_year=2001,
        filter_exists=False,
    ):
        super(AnnualIncidenceAnalyzer, self).__init__(
            working_dir=working_dir,
            filenames=[f"output/MalariaSummaryReport_Yearly_incidence_{start_year}_to_{end_year}.json"]
        )

        self.sweep_variables = sweep_variables or ["Run_Number"]
        self.start_year = start_year
        self.end_year = end_year
        self.filter_exists = filter_exists

    def filter(self, simulation: Simulation):
        if self.filter_exists:
            file = os.path.join(simulation.get_path(), self.filenames[0])
            return os.path.exists(file)
        else:
            return True

    def map(self, data, simulation: Simulation):
        adf = []
        
        fname = self.filenames[0]

        for fname in self.filenames:
            
            age_bins = data[fname]["Metadata"]["Age Bins"]
            years = data[fname]["DataByTime"]["Time Of Report"]
            years = [y / 365 - 1 for y in years]
            
            df = pd.DataFrame.from_dict(data[fname]["DataByTimeAndAgeBins"],orient="columns")[:-1]
            #df["month"] = [[x] * len(age_bins) for x in range(1, len(df) + 1)]
            df["agebin"] = [age_bins] * len(df)
            df["Year"] = [[years[x]] * len(age_bins) for x in range(len(df))]
            df.rename(
               columns={
                   "Annual Clinical Incidence by Age Bin": "Cases",
                   "Average Population by Age Bin": "Pop"
               },
                inplace=True,
            )
            df = df[["agebin", "Year", "Cases", "Pop"]]
            df = df.explode(list(df.columns))
           
            adf = adf + [df]

            adf = pd.concat(adf)

        for sweep_var in self.sweep_variables:
            if sweep_var in simulation.tags.keys():
                adf[sweep_var] = simulation.tags[sweep_var]

        return adf

    def reduce(self, all_data):
        selected = [data for sim, data in all_data.items()]
        if len(selected) == 0:
            print("\nWarning: No data have been returned... Exiting...")
            return

        adf = pd.concat(selected).reset_index(drop=True)
        adf.to_csv(
            (os.path.join(self.working_dir, "ClinicalIncidence_annual.csv")),
            index=False,
        )



class EventReporterAnalyzer(IAnalyzer):
    """
    Pull out the ReportEventRecorder and stack them together.
    """

    def __init__(self, sweep_variables=None, working_dir="./", time_cutoff=0, event_list=["Received_Treatment"], output_filename="events"):
        super(EventReporterAnalyzer, self).__init__(
            working_dir=working_dir, filenames=["output/ReportEventRecorder.csv"]
        )
        self.sweep_variables = sweep_variables
        self.time_cutoff = time_cutoff
        self.event_list = event_list
        self.output_filename = output_filename

    def map(self, data, simulation: Simulation):
        df = data[self.filenames[0]]
        df = df[df["Time"] >= self.time_cutoff].copy()
        df = df[df["Event_Name"].isin(self.event_list)]

        # add tags
        for sweep_var in self.sweep_variables:
            if sweep_var in simulation.tags.keys():
                df[sweep_var] = simulation.tags[sweep_var]

        return df

    def reduce(self, all_data):
        selected = [data for sim, data in all_data.items()]
        if len(selected) == 0:
            print("\nWarning: No data have been returned... Exiting...")
            return

        adf = pd.concat(selected).reset_index(drop=True)
        adf.to_csv(
            (os.path.join(self.working_dir, "".join((self.output_filename,".csv")))),
            index=False,
            index_label=False,
        )


class EventReporterSummaryAnalyzer(IAnalyzer):
    """
    Pull out the ReportEventRecorder and stack them together.
    """

    def __init__(self, sweep_variables=None, working_dir="./", time_cutoff=0, event_list=["Received_Treatment"], output_filename="event_counts"):
        super(EventReporterSummaryAnalyzer, self).__init__(
            working_dir=working_dir, filenames=["output/ReportEventRecorder.csv"]
        )
        self.sweep_variables = sweep_variables
        self.time_cutoff = time_cutoff
        self.event_list = event_list
        self.output_filename = output_filename

    def map(self, data, simulation: Simulation):
        df = data[self.filenames[0]]
        df = df[df["Time"] >= self.time_cutoff].copy()
        df = df[df["Event_Name"].isin(self.event_list)]
        df['Age_Year'] = [age // 365 for age in df['Age']]
        df2 = df.groupby(['Time',"Age_Year",'Event_Name'])['Individual_ID'].agg('count').reset_index()
        df2.rename(columns={"Individual_ID":"Event_Count"})
        # add tags
        for sweep_var in self.sweep_variables:
            if sweep_var in simulation.tags.keys():
                df2[sweep_var] = simulation.tags[sweep_var]

        return df2

    def reduce(self, all_data):
        selected = [data for sim, data in all_data.items()]
        if len(selected) == 0:
            print("\nWarning: No data have been returned... Exiting...")
            return

        adf = pd.concat(selected).reset_index(drop=True)
        adf.to_csv(
            (os.path.join(self.working_dir, "".join((self.output_filename,".csv")))),
            index=False,
            index_label=False,
        )

class NodeDemographicsAnalyzer(IAnalyzer):
    """
    Pull out the NodeDemographicsReport and stack them together.
    """

    def __init__(self, sweep_variables=None, working_dir="./", output_filename="age_population", time_cutoff=0):
        super(NodeDemographicsAnalyzer, self).__init__(
            working_dir=working_dir, filenames=["output/ReportNodeDemographics.csv"]
        )
        self.sweep_variables = sweep_variables
        self.time_cutoff = time_cutoff
        self.output_filename = output_filename

    def map(self, data, simulation: Simulation):
        df = data[self.filenames[0]]
        df = df[df["Time"] >= self.time_cutoff].copy()

        # add tags
        for sweep_var in self.sweep_variables:
            if sweep_var in simulation.tags.keys():
                df[sweep_var] = simulation.tags[sweep_var]

        return df

    def reduce(self, all_data):
        selected = [data for sim, data in all_data.items()]
        if len(selected) == 0:
            print("\nWarning: No data have been returned... Exiting...")
            return

        adf = pd.concat(selected).reset_index(drop=True)
        adf.to_csv(
            (os.path.join(self.working_dir, "".join((self.output_filename,".csv")))),
            index=False,
            index_label=False,
        )

class VectorStatsAnalyzer(IAnalyzer):
    """
    Pull out the Vector Stats reports and stack them together.
    """

    def __init__(self, sweep_variables=None, working_dir="./", start_time=0, end_time=999999):
        super(VectorStatsAnalyzer, self).__init__(
            working_dir=working_dir, filenames=["output/ReportVectorStats.csv"]
        )
        self.sweep_variables = sweep_variables
        self.start_time = start_time
        self.end_time = end_time

    def map(self, data, simulation: Simulation):
        df = data[self.filenames[0]]
        df = df[(df["Time"] >= self.start_time ) & (df['Time'] <= self.end_time)].copy()

        # add tags
        for sweep_var in self.sweep_variables:
            if sweep_var in simulation.tags.keys():
                df[sweep_var] = simulation.tags[sweep_var]

        return df

    def reduce(self, all_data):
        selected = [data for sim, data in all_data.items()]
        if len(selected) == 0:
            print("\nWarning: No data have been returned... Exiting...")
            return

        adf = pd.concat(selected).reset_index(drop=True)
        adf.to_csv(
            (os.path.join(self.working_dir, "VectorStats.csv")),
            index=False,
            index_label=False,
        )


class InsetChartAnalyzer(IAnalyzer):
    """
    Pull out the ReportEventRecorder and stack them together.
    """

    def __init__(
        self, channels, sweep_variables, working_dir="./", start_day=0, end_day=99999
    ):
        super(InsetChartAnalyzer, self).__init__(
            working_dir=working_dir, filenames=["output/InsetChart.json"]
        )
        self.sweep_variables = sweep_variables
        self.channels = channels
        self.start_day = start_day
        self.end_day = end_day

    def map(self, data, simulation: Simulation):
        df = pd.DataFrame(
            {x: data[self.filenames[0]]["Channels"][x]["Data"] for x in self.channels}
        )
        df["time"] = df.index
        df = df[
            (df["time"] >= self.start_day) & (df["time"] <= self.end_day)
        ]
        df["day"] = df["time"] % 365 + 1
        df["year"] = df["time"] // 365

        # add tags
        for sweep_var in self.sweep_variables:
            if sweep_var in simulation.tags.keys():
                df[sweep_var] = simulation.tags[sweep_var]

        return df

    def reduce(self, all_data):
        selected = [data for sim, data in all_data.items()]
        if len(selected) == 0:
            print("\nWarning: No data have been returned... Exiting...")
            return

        adf = pd.concat(selected).reset_index(drop=True)
        adf.to_csv(
            (os.path.join(self.working_dir, "InsetChart.csv")),
            index=False,
            index_label=False,
        )
class PCRAnalyzer(IAnalyzer):
    """
    Pull out the ReportEventRecorder and stack them together.
    """

    def __init__(
        self, channels, sweep_variables, working_dir="./", start_day=0, end_day=99999
    ):
        super(PCRAnalyzer, self).__init__(
            working_dir=working_dir, filenames=["output/InsetChart.json"]
        )
        self.sweep_variables = sweep_variables
        self.channels = channels
        self.start_day = start_day
        self.end_day = end_day

    def map(self, data, simulation: Simulation):
        df = pd.DataFrame(
            {x: data[self.filenames[0]]["Channels"][x]["Data"] for x in self.channels}
        )
        df["time"] = df.index
        df = df[
            (df["time"] >= self.start_day) & (df["time"] <= self.end_day)
        ]
        df["day"] = df["time"] % 365 + 1
        df["year"] = df["time"] // 365

        # add tags
        for sweep_var in self.sweep_variables:
            if sweep_var in simulation.tags.keys():
                df[sweep_var] = simulation.tags[sweep_var]

        return df

    def reduce(self, all_data):
        selected = [data for sim, data in all_data.items()]
        if len(selected) == 0:
            print("\nWarning: No data have been returned... Exiting...")
            return

        adf = pd.concat(selected).reset_index(drop=True)
        adf.to_csv(
            (os.path.join(self.working_dir, "InsetChart_PCR.csv")),
            index=False,
            index_label=False,
        )        
  
class EIRAnalyzer(IAnalyzer):
    """
    Pull out the ReportEventRecorder and stack them together.
    """

    def __init__(
        self, channels, sweep_variables, working_dir="./", start_day=0, end_day=99999
    ):
        super(EIRAnalyzer, self).__init__(
            working_dir=working_dir, filenames=["output/InsetChart.json"]
        )
        self.sweep_variables = sweep_variables
        self.channels = channels
        self.start_day = start_day
        self.end_day = end_day

    def map(self, data, simulation: Simulation):
        df = pd.DataFrame(
            {x: data[self.filenames[0]]["Channels"][x]["Data"] for x in self.channels}
        )
        df["time"] = df.index
        df = df[
            (df["time"] >= self.start_day) & (df["time"] <= self.end_day)
        ]
        df["day"] = df["time"] % 365 + 1
        df["year"] = df["time"] // 365

        # add tags
        for sweep_var in self.sweep_variables:
            if sweep_var in simulation.tags.keys():
                df[sweep_var] = simulation.tags[sweep_var]

        return df

    def reduce(self, all_data):
        selected = [data for sim, data in all_data.items()]
        if len(selected) == 0:
            print("\nWarning: No data have been returned... Exiting...")
            return

        adf = pd.concat(selected).reset_index(drop=True)
        adf.to_csv(
            (os.path.join(self.working_dir, "InsetChart_EIR.csv")),
            index=False,
            index_label=False,
        )
