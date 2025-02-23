import pandas as pd
import os

from wetterdienst.provider.dwd.observation import DwdObservationParameter

from utils import (
    assign_weather_station_to_hospital,
    get_all_weather_data,
    save_git_id,
)

import warnings

warnings.filterwarnings("ignore")


def main(time_interval, save_path, dict_hospital, weather_var_list):
    # load or create hospital_meta
    if os.path.exists(save_path + "/hospital_meta_500missing.csv"):
        hospital_meta = pd.read_csv(
            save_path + "/hospital_meta_500missing.csv", sep=";"
        )

    else:
        # the 500 mising values are only in time interval 2019-01-01 to 2024-01-01
        hospital_meta = assign_weather_station_to_hospital(
            dict_hospital=dict_hospital,
            weather_var_list=weather_var_list,
            time_interval=[
                pd.Timestamp("2019-01-01 00:00:00+00:00"),
                pd.Timestamp("2024-01-01 00:00:00+00:00"),
            ],
        )
        hospital_meta.to_csv(save_path + "/hospital_meta_500missing.csv")

    # create weather data files (full and summary) and save
    get_all_weather_data(
        hospital_meta=hospital_meta,
        variable_list=weather_var_list,
        parth_to_save=save_path,
        time_interval=time_interval,
    )

    save_git_id(path=save_path, file_name="main_weather.py")
    return


if __name__ == "__main__":

    time_interval = [
        pd.Timestamp("2013-01-01 00:00:00+0000"),
        pd.Timestamp("2024-01-01 00:00:00+0000"),
    ]
    save_path = os.getcwd() + "/Weather/data/"

    # All possible Hospitals AKTIN
    dict_hospital = pd.read_csv(
        os.getcwd() + "/Weather/Klinikliste_Adresse_Koordinaten.csv", sep=";"
    )
    dict_hospital.index = dict_hospital["t_klinik"]

    # ------------------------------------------------------------------
    # Weather Data

    # Variable List Weather
    weather_var_list = [
        DwdObservationParameter.HOURLY.TEMPERATURE_AIR_MEAN_200,
        DwdObservationParameter.HOURLY.HUMIDITY,
        DwdObservationParameter.HOURLY.TEMPERATURE_DEW_POINT_MEAN_200,
        DwdObservationParameter.HOURLY.PRESSURE_AIR_SEA_LEVEL,
        DwdObservationParameter.HOURLY.WIND_SPEED,
        DwdObservationParameter.HOURLY.PRECIPITATION_HEIGHT,
    ]

    main(time_interval, save_path, dict_hospital, weather_var_list)
