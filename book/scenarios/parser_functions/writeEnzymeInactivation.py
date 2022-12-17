import pyenzyme as pe
from EnzymePynetics.tools.parameterestimator import ParameterEstimator

def add_enzyme_inactivation_model(
    result: ParameterEstimator, 
    enzmldoc: pe.EnzymeMLDocument, 
    protein_id: str = "p0", 
    reaction_id: str = "r0", 
    mode_name: str = None):
    """Writes enzyme inactivation models to an EnzymeML documante with the respective results from 
    the ParameterEstimator from EnzymePynetics.

    Args:
        result (ParameterEstimator): result from ParameterEstimator
        enzmldoc (pe.EnzymeMLDocument): Enzymeml document from pyEnzyme
        protein_id (str, optional): enzymeml id from the catalyst. Defaults to "p0".
        reaction_id (str, optional): id of the reaction from enzymeml. Defaults to "r0".
        mode_name (str, optional): nome of the model from ParameterEstimator. Defaults to None.
    """

    # define inactive enzyme species
    params = enzmldoc.getProtein('p0').dict()
    inactive_protein = pe.Protein(**params)
    inactive_protein.name += " inactive"
    inactive_protein.init_conc = 0.0
    inactive_protein.constant = False
    inactive_protein.id = "p1"
    enzmldoc.addProtein(inactive_protein)

    # define inactivation reaction
    # Get reaction parameters
    reaction_parameters = enzmldoc.getReaction(reaction_id)
    # Create new reaction
    enzyme_inactivation = pe.EnzymeReaction(
        name="Time-dependent enzyme inactivation",
        reversible=True,
        temperature=reaction_parameters.temperature,
        temperature_unit=reaction_parameters.temperature_unit,
        ph=reaction_parameters.ph)

    # Add reaction species
    enzyme_inactivation.addEduct(protein_id, 1.0, enzmldoc)
    enzyme_inactivation.addProduct("p1", 1.0, enzmldoc)
    enzyme_inactivation_id = enzmldoc.addReaction(enzyme_inactivation)

    # Get k_inact parameter value
    time_unit = result.data.time_unit
    k_inact = pe.enzymeml.models.KineticParameter(
        name="k_inact", 
        value=result.get_model_results(mode_name)["K_ie"].value, 
        unit=f"1 / {time_unit}")

    # Add model result to EnzymeML document
    inactivation_model = pe.KineticModel(
        name="Time-dependent enzyme inactivation",
        equation="- k_inact * p0", 
        parameters=[k_inact],
        enzmldoc=enzmldoc)

    enzyme_inactivation.model = inactivation_model