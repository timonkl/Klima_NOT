import requests
import math
import pandas as pd
import numpy as np
import holidays
import re
import csv
import os
import subprocess

from wetterdienst.provider.dwd.observation import (
    DwdObservationDataset,
    DwdObservationParameter,
    DwdObservationPeriod,
    DwdObservationRequest,
    DwdObservationResolution,
)


def get_distance_in_km(lat1, lat2, long1, long2):
    # Approximation (not accurate for larger distances)
    # see https://www.kompf.de/gps/distcalc.html
    dx = 111.3 * math.cos((lat1 + lat2) / 2 * 0.01745) * (long1 - long2)
    dy = 111.3 * (lat1 - lat2)
    return math.sqrt(dx**2 + dy**2)


def get_weather_station_meta_data():
    weather_data = DwdObservationRequest(
        parameter=DwdObservationDataset.TEMPERATURE_AIR,
        resolution=DwdObservationResolution.MINUTE_10,
        period=DwdObservationPeriod.HISTORICAL,
    ).all()
    weather_data.df.head()
    return weather_data.df.to_pandas()


def find_nearest_weather_station(lat, long, list_station, weather_station_meta):
    dist = [
        get_distance_in_km(
            lat1=lat, lat2=row["latitude"], long1=long, long2=row["longitude"]
        )
        for idx, row in list_station.iterrows()
    ]

    sorted_idx_distance = sorted(range(len(dist)), key=lambda k: dist[k])
    sorted_dist = np.array(dist)[sorted_idx_distance]

    # get the Weather statioin IDs to the corresponding INDEX
    id_sorted_distance = [
        int(weather_station_meta["station_id"][sorted_idx_distance[i]])
        for i in range(len(sorted_idx_distance))
    ]

    df = pd.DataFrame(
        {
            "station_id": id_sorted_distance,
            "station_idx": sorted_idx_distance,
            "distance_to_hospital": sorted_dist,
        }
    )
    return df


def assign_weather_station_to_hospital(
    dict_hospital,
    weather_var_list,
    time_interval=[
        pd.Timestamp("2019-01-01 00:00:00+00:00"),
        pd.Timestamp("2024-01-01 00:00:00+00:00"),
    ],
    max_missing_values=500,
    max_dist_weather_station=50,
):

    # get the meta data from the hospitals
    hospital_meta = dict_hospital.copy()

    # add empty Columns for meta data of the assigned weather stations
    hospital_meta["weather_station_id"] = ""
    hospital_meta["weather_station_latitude"] = ""
    hospital_meta["weather_station_longitude"] = ""
    hospital_meta["weather_station_height"] = ""
    hospital_meta["weather_station_distance_hospital"] = ""
    hospital_meta["weather_station_start_date"] = ""
    hospital_meta["weather_station_end_date"] = ""
    hospital_meta["weather_station_missing_values"] = ""
    hospital_meta["weather_station_comment"] = ""

    # get meta data of all weather stations
    weather_station_meta = get_weather_station_meta_data()

    # iterate over all hospitals
    for hospital_id, hospital in hospital_meta["klinik_name"].items():

        print("Processing Hospital: " + str(hospital_meta["klinik_name"][hospital_id]))

        weather_station_sorted_by_distance = find_nearest_weather_station(
            lat=hospital_meta["Latitude"][hospital_id],
            long=hospital_meta["Longitude"][hospital_id],
            list_station=weather_station_meta,
            weather_station_meta=weather_station_meta,
        )

        for i in range(weather_station_sorted_by_distance.shape[0]):

            # to get ride of the long names
            weather_id = weather_station_sorted_by_distance["station_id"][i]
            weather_idx = weather_station_sorted_by_distance["station_idx"][i]
            distance = weather_station_sorted_by_distance["distance_to_hospital"][i]

            if distance > max_dist_weather_station:
                print(
                    f"There is no weather station for {hospital} in {max_dist_weather_station} km distance"
                )
                break

            table = get_weather_data(weather_var_list, station_id=weather_id)

            # check if the weather station is available
            if table is None:
                comment = (
                    hospital_meta["weather_station_comment"][hospital_id]
                    + str(weather_id)
                    + ", "
                )
                hospital_meta.loc[hospital_id, "weather_station_comment"] = comment
                continue

            # Find the missing dates
            full_date_range = pd.date_range(
                start=time_interval[0], end=time_interval[1], freq="h"
            )

            missing_dates = full_date_range.difference(table["date"])

            # Skip to the next station if missing values exceed the limit
            if missing_dates.to_series().count() > max_missing_values:
                comment = (
                    hospital_meta["weather_station_comment"][hospital_id]
                    + str(weather_id)
                    + ", "
                )
                hospital_meta.loc[hospital_id, "weather_station_comment"] = comment
                continue

            # add stats of missing values
            missing_dates = full_date_range.difference(table["date"])
            missing_dates_per_year = (
                missing_dates.to_series().groupby(missing_dates.year).count()
            )
            missing_per_year = ""
            for year, count_missing in missing_dates_per_year.items():
                missing_per_year += " " + str(year) + ": " + str(count_missing) + ","

            # add data to the table
            hospital_meta.loc[hospital_id, "weather_station_id"] = weather_station_meta[
                "station_id"
            ][weather_idx]
            hospital_meta.loc[hospital_id, "weather_station_latitude"] = (
                weather_station_meta["latitude"][weather_idx]
            )
            hospital_meta.loc[hospital_id, "weather_station_longitude"] = (
                weather_station_meta["longitude"][weather_idx]
            )
            hospital_meta.loc[hospital_id, "weather_station_height"] = (
                weather_station_meta["height"][weather_idx]
            )
            hospital_meta.loc[hospital_id, "weather_station_distance_hospital"] = (
                distance
            )
            hospital_meta.loc[hospital_id, "weather_station_start_date"] = table[
                "date"
            ].min()
            hospital_meta.loc[hospital_id, "weather_station_end_date"] = table[
                "date"
            ].max()
            hospital_meta.loc[hospital_id, "weather_station_missing_values"] = (
                missing_per_year
            )

            break

    return hospital_meta


def get_weather_data(
    variable_list,
    station_id=1048,
    time_interval=[
        pd.Timestamp("2019-01-01 00:00:00+00:00"),
        pd.Timestamp("2024-01-01 00:00:00+00:00"),
    ],
    max_iterations=10,
):
    for itr, itr_var in enumerate(variable_list):

        if itr == 0:
            table_0 = (
                DwdObservationRequest(
                    parameter=itr_var,
                    resolution=DwdObservationResolution.HOURLY,
                    period=DwdObservationPeriod.HISTORICAL,
                )
                .filter_by_station_id(station_id=station_id)
                .values.all()
            ).df.to_pandas()

            if table_0.empty:
                print("There are no Data for Weather station {}".format(station_id))
                return

            variable_name = table_0["parameter"][0]
            table_0 = table_0[["station_id", "date", "value"]]
            table_0.rename(columns={"value": variable_name}, inplace=True)

        else:
            table_1 = (
                DwdObservationRequest(
                    parameter=itr_var,
                    resolution=DwdObservationResolution.HOURLY,
                    period=DwdObservationPeriod.HISTORICAL,
                )
                .filter_by_station_id(station_id=station_id)
                .values.all()
            ).df.to_pandas()

            if table_1.empty:
                print("There are no Data for Weather station {}".format(station_id))
                return

            variable_name = table_1["parameter"][0]
            table_1 = table_1[["station_id", "date", "value"]]
            table_1.rename(columns={"value": variable_name}, inplace=True)
            table_0 = pd.merge(table_0, table_1)

    table_0["date"] = pd.to_datetime(table_0["date"])

    if time_interval:
        mask = (table_0["date"] >= time_interval[0]) & (
            table_0["date"] <= time_interval[1]
        )
        table_0 = table_0.loc[mask].reset_index(drop=True)

    # from Kelvin to Celsius
    if "temperature_air_mean_200" in table_0.keys():
        table_0["temperature_air_mean_200"] = (
            table_0["temperature_air_mean_200"] - 273.15
        )

    if "temperature_dew_point_mean_200" in table_0.keys():
        table_0["temperature_dew_point_mean_200"] = (
            table_0["temperature_dew_point_mean_200"] - 273.15
        )

    return table_0


def get_all_weather_data(
    hospital_meta,
    variable_list,
    parth_to_save,
    time_interval=[
        pd.Timestamp("2019-01-01 00:00:00+00:00"),
        pd.Timestamp("2024-01-01 00:00:00+00:00"),
    ],
    max_missing_entries=5,
):

    bundeslaender = [
        "Brandenburg",
        "Berlin",
        "Baden-Württemberg",
        "Bayern",
        "Bremen",
        "Hessen",
        "Hamburg",
        "Mecklenburg-Vorpommern",
        "Niedersachsen",
        "Nordrhein-Westfalen",
        "Rheinland-Pfalz",
        "Schleswig-Holstein",
        "Saarland",
        "Sachsen",
        "Sachsen-Anhalt",
        "Thüringen",
    ]
    bundeslaender_codes = [
        "BB",
        "BE",
        "BW",
        "BY",
        "HB",
        "HE",
        "HH",
        "MV",
        "NI",
        "NW",
        "RP",
        "SH",
        "SL",
        "SN",
        "ST",
        "TH",
    ]

    for idx_hospital in range(hospital_meta.shape[0]):

        """
        # nice for debug
        if idx_hospital > 2:
            break
        """
        table = get_weather_data(
            variable_list=variable_list,
            station_id=hospital_meta.loc[idx_hospital, "weather_station_id"],
            time_interval=time_interval,
        )

        table["t_klinik"] = hospital_meta.loc[idx_hospital, "t_klinik"]

        # drop redundant values
        if table["station_id"].nunique() == 1:
            # table = table.drop(columns=["station_id"])
            pass
        else:
            raise ValueError("Station id changed over time!")

        # add weekday
        table.insert(1, "weekday", table["date"].dt.weekday)  # Monday=0, Sunday=6

        # add holiday
        bundesland = hospital_meta.loc[idx_hospital, "Bundesland"]
        de_holidays = holidays.Germany(
            state=bundeslaender_codes[
                np.where(np.array(bundeslaender) == bundesland)[0][0]
            ]
        )
        table["holiday"] = table["date"].dt.date.apply(lambda _: _ in de_holidays)

        # save full table
        table.reset_index()
        table.to_csv(
            parth_to_save
            + "/hospital_id_"
            + str(hospital_meta.loc[idx_hospital, "t_klinik"])
            + ".csv"
        )

        #########
        # Summary

        # temperature
        nan_counts = table.groupby(table["date"].dt.date)[
            "temperature_air_mean_200"
        ].apply(lambda x: x.isna().sum())
        days_to_remove = nan_counts[nan_counts > max_missing_entries].index
        table_temperature = (
            table.groupby(table["date"].dt.date)["temperature_air_mean_200"]
            .agg(["mean", "max", "min"])
            .reset_index()
        )
        table_temperature["max_difference"] = (
            table_temperature["max"] - table_temperature["min"]
        )
        table_temperature.columns = [
            "date",
            "temperature_air_mean",
            "temperature_air_max",
            "temperature_air_min",
            "temperature_air_temp_max_difference",
        ]
        for day in days_to_remove:
            table_temperature.loc[
                table_temperature["date"] == day,
                [
                    "temperature_air_mean",
                    "temperature_air_max",
                    "temperature_air_min",
                    "temperature_air_temp_max_difference",
                ],
            ] = np.nan

        # humidity
        nan_counts = table.groupby(table["date"].dt.date)["humidity"].apply(
            lambda x: x.isna().sum()
        )
        days_to_remove = nan_counts[nan_counts > max_missing_entries].index
        table_humidity = (
            table.groupby(table["date"].dt.date)["humidity"]
            .agg(["mean", "max", "min"])
            .reset_index()
        )
        table_humidity.columns = [
            "date",
            "humidity_mean",
            "humidity_max",
            "humidity_min",
        ]
        for day in days_to_remove:
            table_humidity.loc[
                table_humidity["date"] == day,
                [
                    "humidity_mean",
                    "humidity_max",
                    "humidity_min",
                ],
            ] = np.nan

        # dew_point
        nan_counts = table.groupby(table["date"].dt.date)[
            "temperature_dew_point_mean_200"
        ].apply(lambda x: x.isna().sum())
        days_to_remove = nan_counts[nan_counts > max_missing_entries].index
        table_dew_point = (
            table.groupby(table["date"].dt.date)["temperature_dew_point_mean_200"]
            .agg(["mean", "max", "min"])
            .reset_index()
        )
        table_dew_point.columns = [
            "date",
            "temperature_dew_point_mean",
            "temperature_dew_point_max",
            "temperature_dew_point_min",
        ]
        for day in days_to_remove:
            table_dew_point.loc[
                table_dew_point["date"] == day,
                [
                    "temperature_dew_point_mean",
                    "temperature_dew_point_max",
                    "temperature_dew_point_min",
                ],
            ] = np.nan

        # pressure_air
        nan_counts = table.groupby(table["date"].dt.date)[
            "pressure_air_sea_level"
        ].apply(lambda x: x.isna().sum())
        days_to_remove = nan_counts[nan_counts > max_missing_entries].index
        table_pressure_air = (
            table.groupby(table["date"].dt.date)["pressure_air_sea_level"]
            .agg(["mean", "max", "min"])
            .reset_index()
        )
        table_pressure_air.columns = [
            "date",
            "pressure_air_sea_level_mean",
            "pressure_air_sea_level_max",
            "pressure_air_sea_level_min",
        ]
        for day in days_to_remove:
            table_pressure_air.loc[
                table_dew_point["date"] == day,
                [
                    "pressure_air_sea_level_mean",
                    "pressure_air_sea_level_max",
                    "pressure_air_sea_level_min",
                ],
            ] = np.nan

        # wind_speed
        nan_counts = table.groupby(table["date"].dt.date)["wind_speed"].apply(
            lambda x: x.isna().sum()
        )
        days_to_remove = nan_counts[nan_counts > max_missing_entries].index
        table_wind_speed = (
            table.groupby(table["date"].dt.date)["wind_speed"]
            .agg(["mean", "max", "min"])
            .reset_index()
        )
        table_wind_speed.columns = [
            "date",
            "wind_speed_mean",
            "wind_speed_max",
            "wind_speed_min",
        ]
        for day in days_to_remove:
            table_wind_speed.loc[
                table_dew_point["date"] == day,
                [
                    "wind_speed_mean",
                    "wind_speed_max",
                    "wind_speed_min",
                ],
            ] = np.nan

        # precipitation_height
        nan_counts = table.groupby(table["date"].dt.date)["precipitation_height"].apply(
            lambda x: x.isna().sum()
        )
        days_to_remove = nan_counts[nan_counts > max_missing_entries].index
        table_precipitation = (
            table.groupby(table["date"].dt.date)["precipitation_height"]
            .sum()
            .reset_index()
        )
        table_precipitation.columns = ["date", "precipitation_height_sum"]
        for day in days_to_remove:
            table_precipitation.loc[
                table_dew_point["date"] == day,
                [
                    "precipitation_height_sum",
                ],
            ] = np.nan

        table_summary = pd.merge(
            table_temperature, table_humidity, on="date", how="outer"
        )
        table_summary = pd.merge(table_summary, table_dew_point, on="date", how="outer")
        table_summary = pd.merge(
            table_summary, table_pressure_air, on="date", how="outer"
        )
        table_summary = pd.merge(
            table_summary, table_wind_speed, on="date", how="outer"
        )
        table_summary = pd.merge(
            table_summary, table_precipitation, on="date", how="outer"
        )

        table_summary["date"] = pd.to_datetime(table_summary["date"])
        # add weekday
        table_summary.insert(
            1, "weekday", table_summary["date"].dt.weekday
        )  # Monday=0, Sunday=6

        # add holiday
        bundesland = hospital_meta.loc[idx_hospital, "Bundesland"]
        de_holidays = holidays.Germany(
            state=bundeslaender_codes[
                np.where(np.array(bundeslaender) == bundesland)[0][0]
            ]
        )
        table_summary["holiday"] = table_summary["date"].dt.date.apply(
            lambda _: _ in de_holidays
        )

        table_summary.reset_index()

        table_summary.to_csv(
            parth_to_save
            + "/summary_hospital_id_"
            + str(hospital_meta.loc[idx_hospital, "t_klinik"])
            + ".csv"
        )


# Air Quality Data Umwelt Bundesamt
# ----------------------------------------------------------------------------------------------------------------------


def get_air_data(
    date_from="2023-03-03",
    date_to="2023-03-04",
    station_id="DEST103",
    component="1",
    scope="1",
):
    """
    get the Air Data from Umwelt Bundesamt, description see attached PDF
    """

    component_dict = {
        "1": ["Feinstaub", "PM10"],
        "2": ["Kohlenmonoxid", "CO"],
        "3": ["Ozon", "O3"],
        "4": ["Schwefeldioxid", "SO2"],
        "5": ["Stickstoffdioxid", "NO2"],
        "6": ["Blei im Feinstaub", "PB"],
        "7": ["Benzo(a)pyrene im Feinstaub", "BaP"],
        "8": ["Benzol", "C6H6"],
        "9": ["Feinstaub", "PM25"],
        "10": ["Arsen im Feinstaub", "As"],
        "11": ["Cadmium im Feinstaub", "Cd"],
        "12": ["Nickel im Feinstaub", "Ni"],
    }

    scope_dict = {
        "1": ["Tagesmittel", "1TMW"],
        "2": ["Ein-Stunden-Mittelwert", "1SMW"],
        "3": ["Ein-Stunden-Tagesmaxima", "1SMW_MAX"],
        "4": ["Acht-Stunden-Mittelwert", "8SMW"],
        "5": ["Acht-Stunden-Tagesmaxima", "8SMW_MAX"],
        "6": ["stündlich gleitendes Tagesmittel", "1TMWGL"],
    }

    # all stations
    # response = requests.get('https://www.umweltbundesamt.de/api/air_data/v3/meta/json?use=airquality&' +
    #                         'date_from=2023-01-01&date_to=2023-06-30&time_from=1&time_to=24&lang=de')
    # for key in response['stations'].keys():
    #     print(response['stations'][key][2])

    link_to_api = (
        "https://www.umweltbundesamt.de/api/air_data/v3/measures/json?date_from="
        + date_from
        + "&date_to="
        + date_to
        + "&time_from=1&time_to=24&station="
        + station_id
        + "&component="
        + component
        + "&scope="
        + scope
    )

    response = requests.get(link_to_api)
    return response.json()


def get_air_meta(
    date_from="2023-03-03",
    date_to="2023-03-04",
    component=[],
    scope="1",
    save_path=False,
):

    link_to_api = (
        "https://www.umweltbundesamt.de/api/air_data/v3/meta/json?use=airquality&date_from="
        + date_from
        + "&date_to="
        + date_to
        + "&time_from=1&time_to=24"
    )

    response = requests.get(link_to_api)

    air_meta = pd.DataFrame.from_dict(response.json()["stations"])
    air_meta = air_meta.T

    air_meta.columns = [
        "station_number",
        "station_id",
        "station_name",
        "city",
        " ",
        "start_date",
        "end_date",
        "longitude",
        "latitude",
        " ",
        " ",
        " ",
        "state_code",
        "state",
        "area_type",
        "zone",
        "influence",
        "street",
        "house_number",
        "postal_code",
    ]

    air_meta["latitude"] = air_meta["latitude"].astype(float)
    air_meta["longitude"] = air_meta["longitude"].astype(float)

    # check the provided components
    if len(component) > 0:
        air_meta["component"] = ""
        for idx, row in air_meta.iterrows():

            print(row["station_id"])
            component_provided = ""
            for itr_component in component:
                try:
                    response = get_air_data(
                        date_from=date_from,
                        date_to=date_to,
                        station_id=row["station_id"],
                        component=itr_component,
                        scope=scope,
                    )

                    if response["data"]:
                        component_provided += str(itr_component) + ", "
                except:
                    component_provided = "ERROR"
                    break

            air_meta.loc[idx, "component"] = component_provided

            if save_path:
                air_meta.to_csv(save_path)

    return air_meta


def find_nearest_air_station(lat, long, dict_station):
    dist = [
        get_distance_in_km(
            lat1=lat, lat2=row["latitude"], long1=long, long2=row["longitude"]
        )
        for idx, row in dict_station.iterrows()
    ]

    sorted_idx_distance = sorted(range(len(dist)), key=lambda k: dist[k])
    sorted_dist = np.array(dist)[sorted_idx_distance]

    id_sorted = [
        str(dict_station["station_id"][sorted_idx_distance[i]])
        for i in range(len(sorted_idx_distance))
    ]

    component_sorted = [
        str(dict_station["component"][sorted_idx_distance[i]])
        for i in range(len(sorted_idx_distance))
    ]

    df = pd.DataFrame(
        {
            "station_id": id_sorted,
            "distance_to_hospital": sorted_dist,
            "component": component_sorted,
        }
    )
    return df


def get_all_air_data(
    save_path,
    hospital_air_meta,
    air_meta,
    component,
    scope,
    max_distance,
    date_from,
    date_to,
    max_missing_entries=5,
):

    component_dict = {
        "1": ["Feinstaub", "PM10"],
        "2": ["Kohlenmonoxid", "CO"],
        "3": ["Ozon", "O3"],
        "4": ["Schwefeldioxid", "SO2"],
        "5": ["Stickstoffdioxid", "NO2"],
        "6": ["Blei im Feinstaub", "PB"],
        "7": ["Benzo(a)pyrene im Feinstaub", "BaP"],
        "8": ["Benzol", "C6H6"],
        "9": ["Feinstaub", "PM25"],
        "10": ["Arsen im Feinstaub", "As"],
        "11": ["Cadmium im Feinstaub", "Cd"],
        "12": ["Nickel im Feinstaub", "Ni"],
    }

    if os.path.exists(save_path + "meta_air" + ".csv"):
        hospital_air_meta = pd.read_csv(save_path + "meta_air" + ".csv", sep=",")
        print("meta_air already exists")
    else:
        hospital_air_meta["air_id"] = 0
        hospital_air_meta["air_latitude"] = ""
        hospital_air_meta["air_longitude"] = ""
        hospital_air_meta["air_distance"] = ""
        hospital_air_meta["air_component"] = ""

    for idx_hospital, row_hospital in hospital_air_meta.iterrows():

        # check if air data for hospital already exists
        if isinstance(row_hospital["air_id"], str) and row_hospital["air_id"] != "0":
            print(f"Skip {row_hospital["air_id"]}")
            continue

        dist_air = find_nearest_air_station(
            lat=row_hospital["Latitude"],
            long=row_hospital["Longitude"],
            dict_station=air_meta,
        )

        for idx, row in dist_air.iterrows():

            # Break if there is no air data in max_distance
            if row["distance_to_hospital"] > max_distance:
                print(
                    f"No air station in distance to Hospital {row_hospital['t_klinik']}"
                )
                break

            if all(element in row["component"] for element in component):

                merged_data = pd.DataFrame(columns=["Date"])
                for i, component_itr in enumerate(component):

                    try:
                        air_data = get_air_data(
                            date_from=date_from,
                            date_to=date_to,
                            station_id=row["station_id"],
                            component=component_itr,
                            scope=scope,
                        )

                        data = []
                        for date, values in air_data["data"][
                            list(air_data["data"].keys())[0]
                        ].items():
                            data.append({"Date": date, "Value": values[2]})

                        # Create DataFrame
                        df = pd.DataFrame(data)

                        df = df.rename(
                            columns={"Value": component_dict[component_itr][1]}
                        )

                        merged_data = pd.merge(merged_data, df, on="Date", how="outer")

                        # Update meta
                        hospital_air_meta.at[idx_hospital, "air_id"] = row["station_id"]
                        hospital_air_meta.at[idx_hospital, "air_latitude"] = (
                            air_meta.loc[
                                air_meta["station_id"] == row["station_id"], "latitude"
                            ].iloc[0]
                        )

                        hospital_air_meta.at[idx_hospital, "air_longitude"] = (
                            air_meta.loc[
                                air_meta["station_id"] == row["station_id"], "longitude"
                            ].iloc[0]
                        )
                        hospital_air_meta.at[idx_hospital, "air_distance"] = row[
                            "distance_to_hospital"
                        ]
                        hospital_air_meta.at[idx_hospital, "air_component"] = row[
                            "component"
                        ]

                    except:
                        # change to hospital id!
                        print(
                            f"No data for component {component_itr} in station {row['station_id']}"
                        )

                # save data
                merged_data.to_csv(
                    save_path + "hospital_id_" + str(row_hospital["t_klinik"]) + ".csv"
                )

                merged_data["Date"] = pd.to_datetime(merged_data["Date"])

                # set all values for days with more than max_missing_entries in the current column to NaN
                columns = merged_data.columns.difference(["Date"])
                for col in columns:
                    nan_counts = merged_data.groupby(merged_data["Date"].dt.date)[
                        col
                    ].apply(lambda x: x.isna().sum())
                    days_to_remove = nan_counts[nan_counts > max_missing_entries].index

                    for day in days_to_remove:
                        merged_data.loc[
                            merged_data["Date"].dt.date == day,
                            col,
                        ] = np.nan

                merged_data["Date"] = pd.to_datetime(merged_data["Date"])
                merged_data.set_index("Date", inplace=True)
                # Resample by day and calculate mean, min, max
                daily_stats = merged_data.resample("D").agg(["mean", "min", "max"])

                # Reset column names for better readability
                daily_stats.columns = [
                    "_".join(col).strip() for col in daily_stats.columns
                ]

                daily_stats.to_csv(
                    save_path
                    + "summary_hospital_id_"
                    + str(row_hospital["t_klinik"])
                    + ".csv"
                )

                break

    # save meta
    hospital_air_meta.to_csv(save_path + "meta_air" + ".csv")

    return


# GitHub ID
# ----------------------------------------------------------------------------------------------------------------------


def save_git_id(path, file_name):

    try:
        commit_id = (
            subprocess.check_output(["git", "rev-parse", "HEAD"])
            .strip()
            .decode("utf-8")
        )

        with open(path + "git_report.txt", "a") as file:
            file.write(f"Latest Git Commit ID for {file_name}: {commit_id}\n")

    except subprocess.CalledProcessError:
        print("Error: Not a valid Git repository or Git is not installed.")
