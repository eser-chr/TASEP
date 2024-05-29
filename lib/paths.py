from dataclasses import dataclass

@dataclass
class Paths:
    data_folder:str = "./DATA"
    last_id_path = "./id"
    catalog = "./catalog.csv"
    library = "./LIBS"