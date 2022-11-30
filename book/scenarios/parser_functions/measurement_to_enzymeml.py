from typing import List
import pyenzyme as pe

from .utility import read_photometer, tab_split, string_to_float, to_seconds


def measurement_data_to_EnzymeML(
    template_path: str,
    measurement_data: dict,
    species_ids: List[str],
    data_unit: str,
    time_unit: str
    ) -> pe.EnzymeMLDocument:

    enzmldoc: pe.EnzymeMLDocument = pe.EnzymeMLDocument.fromTemplate(template_path)
    enzmldoc.name = measurement_data["name"]
    enzmldoc.created = measurement_data["date"]
    enzmldoc.reaction_dict["r0"].ph = measurement_data["pH"]
    enzmldoc.reaction_dict["r0"].temperature = measurement_data["temperature"]
    enzmldoc.reaction_dict["r0"].temperature_unit = "C"
    
    for species in species_ids:

        for IDs, concentration in zip(enzmldoc.measurement_dict.keys(), measurement_data[species]):
            for counter, replicate in enumerate(concentration["data"]):
                
                enzmldoc.getMeasurement(IDs).addReplicates(pe.Replicate(
                    id=f"Measurement{counter}",
                    species_id=species,
                    data=replicate["replicate"],
                    data_unit=data_unit,
                    time=measurement_data["time"],
                    time_unit=time_unit), enzmldoc)

    return enzmldoc