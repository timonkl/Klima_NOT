import pandas as pd
import os

from utils import get_all_air_data, get_air_meta, save_git_id


def main(save_path, hospital_meta, component, scope, max_distance, date_from, date_to):
    if os.path.exists(save_path + "meta_air_scope_" + scope + ".csv"):
        air_meta = pd.read_csv(save_path + "meta_air_scope_" + scope + ".csv", sep=",")
    else:
        air_meta = get_air_meta(
            component=[str(i) for i in range(1, 13)],
            scope=scope,
            save_path=save_path + "meta_air_scope_" + scope + ".csv",
        )

    get_all_air_data(
        save_path=save_path,
        hospital_air_meta=hospital_meta,
        air_meta=air_meta,
        component=component,
        scope=scope,
        max_distance=max_distance,
        date_from=date_from,
        date_to=date_to,
    )

    save_git_id(path=save_path, file_name="main_air.py")

if __name__ == "__main__":

    save_path = os.getcwd() + "/Weather/air/"
    component = ["1", "9", "3", "5"]
    scope = "2"
    max_distance = 30
    date_from = "2019-01-01"
    date_to = "2024-01-01"

    # All possible Hospitals AKTIN
    dict_hospital = pd.read_csv(
        os.getcwd() + "/Weather/Klinikliste_Adresse_Koordinaten.csv", sep=";"
    )
    dict_hospital.index = dict_hospital["t_klinik"]

    main(
        save_path=save_path,
        component=component,
        scope=scope,
        max_distance=max_distance,
        date_from=date_from,
        date_to=date_to,
        hospital_meta=dict_hospital,
    )
