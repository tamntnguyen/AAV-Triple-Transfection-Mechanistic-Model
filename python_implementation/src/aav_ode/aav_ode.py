#!/usr/bin/env python

"""
Mechanistic model using ordinary differential equations (ODE)
to predict the dynamics of adenovirus-associated viral
replication and production in vitro.

Nguyen TNT, Sha S, Hong MS, Maloney AJ, Barone PW, Neufeld C, 
Wolfrum J, Springs SL, Sinskey AJ, Braatz RD. 
Mechanistic model for production of recombinant adeno-associated 
virus via triple transfection of HEK293 cells. 
Mol Ther Methods Clin Dev. 2021 Apr 16;21:642-655. 
doi: 10.1016/j.omt m.2021.04.006. 
PMID: 34095346; 
PMCID: PMC8143981.
https://www.cell.com/molecular-therapy-family/methods/fulltext/S2329-0501(21)00072-3

Original MATLAB source code: 
https://github.com/tamntnguyen/AAV-Triple-Transfection-Mechanistic-Model

"""

from enum import IntEnum, unique

import numpy as np  # pip install numpy
import click  # pip install click

import plotly.express as px  # pip install plotly
import plotly.graph_objects as go
from plotly.subplots import make_subplots

# For saving images you'll also need Kaleido
# pip install kaleido

from scipy.integrate import solve_ivp  # pip install scipy

PACKAGE_NAME = "aav_ode"


@unique
class kType(IntEnum):
    """Constants for the model

    ' Figure 2.
    Helper plasmid (pHelper) activates expression of the rep/
    cap gene on the packaging plasmid (pPackaging) and
    synthesis (kCap_syn and kRep_syn) of viral protein (VP)
    and Rep protein (RepProtein). Capsid proteins are
    assembled (kassembly) into empty capsids in the nucleus
    (EmptyCapNuc), and each capsid particle consists of 60
    protein subunits. With helper functions from the helper
    plasmid, the Rep protein replicates (kDNA_rep) viral DNA
    (vDNA) from the transgene vector plasmid (pVector). Rep
    proteins dock on empty capsids (kRep_bind_capsid) to form
    intermediate complexes (CapRepcomplex), which then
    interact with viral DNA and encapsidate them inside capsids
    (kDNA_pack) at a 1:1 ratio to form full virions inside the nucleus
    (RepRCcomplex). Regardless of their contents (empty or
    full), capsids can be secreted out of the nucleus (ksecrete) into
    the cytosol (EmptyCapCyto and FullCapCyto). Rep protein
    binding (kRep_bind_plasmid) forming bounded plasmids
    (RepRCcomplex) negatively regulates expression of the
    rep/cap gene on the packaging plasmid. Degradation of
    proteins (kRep_protein_degrade and kVP_degrade) is possible
    during the process. '

    """

    kUptake = 0
    kEscape = 1
    kNuclearEntry = 2
    kPlasmidDeg = 3
    kRepSyn = 4
    kCapSyn = 5
    kProteinDegRep = 6  # K_Rep_protein_degrade
    kProteinDegCap = 7  # kVP_degrade
    kDNArep = 8
    kBindRCplasmid = 9
    kAssembly = 10
    kBindCapsid = 11
    kPack = 12  # k_DNA_pack
    kSecrete = 13


@unique
class outputTypes(IntEnum):
    """Types of output from the model"""

    pRC_extracellular = 0
    pRC_endosomal = 1
    pRC_cytosol = 2
    Packaging_Plasmid = 3  # aka "pRC nucleus"
    pVector_extracellular = 4
    pVector_endosomal = 5
    pVector_cytosol = 6
    Vector_Plasmid = 7  # aka "pVector nucleus"
    pHelper_extracellular = 8
    pHelper_endosomal = 9
    pHelper_cytosol = 10
    Helper_Plasmid = 11  # aka "pHelper nucleus"
    Rep_Protein = 12
    VP = 13  # aka "cap protein"
    vDNA = 14
    Empty_Capsid_nucleus = 15
    Full_Capsid_nucleus = 16
    Empty_Capsid_cytosol = 17
    Full_Capsid_cytosol = 18
    Rep_bound_Packaging_Plasmid = 19
    Rep_bound_Empty_Plasmid = 20


class MEASURED_DATA:
    """
    Experimental data from the original paper
    https://www.cell.com/molecular-therapy-family/methods/fulltext/S2329-0501(21)00072-3

    ' Model simulation with estimated parameters provides a good fit with
    experimental data (Figure 4). Interestingly, the plasmid kuptake
    estimated from in-house experimental data is significantly lower
    than that estimated from literature data. In-house experimental data
    show that viral DNA started replicating around 12 hours
    post transfection (hpt) and doubled every ∼12 h (Figure 4A).
    The total amount of capsid produced per cell increased sharply
    at first, with ∼80% of total capsids produced within the
    first 24 h and then plateaued. Total capsid particles and full
    virion production display a similar trend and saturate at
    a slightly later time point between 24 and 36 h; this trend is
    consistent with data from another transfection study, where the
    majority of full virions was produced in the first 24 h.
    Given the assumption that newly synthesized capsid proteins
    are rapidly assembled, the plateau in capsid particle production
    would be caused by saturation in VP production regulated by total
    Rep protein production. The amount of full virion only
    accounts for a consistent 2%–3% of the total capsid produced
    (Figure 4B). Although most of viral DNA was encapsidated at
    24 hpt, and DNA replication continued at an exponential rate,
    the rate of encapsidated DNA did not follow DNA replication
    and decreased significantly beyond that time point.
    By 48 hpt, only 26% of total viral DNA made it inside capsid
    particles (Figure 4A). '

    """

    MEDIA_EXCHANGE_TIME = 6.0  # Hours from transfection to media exchange
    TOTAL_TIME = 60.0  # Hours of experiment after transfection

    # TODO: Is this [rep/cap plasmid, GOI plasmid, helper plasmid]?
    INITIAL_CONCENTRATIONS = [7.6e4, 7.6e4, 7.6e4]  # Initial concentrations

    kUptake = 1.1888e-03  # per hour
    kEscape = 6.0000e-01  # per hour
    kNuclearEntry = 4.3000e-03  # per hour
    kPlasmidDeg = 1.9500e-02  # per hour
    kRepSyn = 6.4891e04  # molecule per cell per hour
    kCapSyn = 6.4891e04  # molecule per cell per hour
    kProteinDegRep = 2.4500e-02  # per hour
    kProteinDegCap = 2.7730e-01  # per hour
    kDNArep = 3.0990e-07  # cell^2 per molecule^2 per hour
    kBindRCplasmid = 3.3829e-05  # molecule per cell per hour
    kAssembly = 1.0000e05  # per hour
    kBindCapsid = 5.6979e-03  # molecule per cell per hour
    kPack = 8.1784e-02  # moelcule per cell per hour
    kSecrete = 1.0000e05  # per hour

    data_Total_DNA = np.array(
        [
            [10, 5.6604e08, 5.8795e08, 4.9559e08],
            [24, 9.8793e08, 1.3904e09, 1.0795e09],
            [34, 2.2493e09, 2.221e09, 1.791e09],
            [46, 4.5122e09, 4.5409e09, 5.3871e09],
        ]
    )

    data_full_cap = np.array(
        [
            [24, 9.6824e08, 8.318e08, 7.8574e08],
            [34, 1.1707e09, 1.0381e09, 9.6213e08],
            [46, 1.2472e09, 1.0446e09, 9.6824e08],
        ]
    )

    data_cap = np.array(
        [
            [10, 5.85e08],
            [10, 5.84e08],
            [24, 3.62e10],
            [24, 4.03e10],
            [34, 4.11e10],
            [34, 3.59e10],
            [46, 5.58e10],
            [46, 5.29e10],
            [55, 4.52e10],
            [55, 4.77e10],
        ]
    )


class Experiments:
    """Experiments class object

    Models the adeno-associated virus (AAV) transfection.
    Based on the paper:

    Nguyen TNT, Sha S, Hong MS, Maloney AJ, Barone PW, Neufeld C,
    Wolfrum J, Springs SL, Sinskey AJ, Braatz RD.
    Mechanistic model for production of recombinant adeno-associated
    virus via triple transfection of HEK293 cells.
    Mol Ther Methods Clin Dev. 2021 Apr 16;21:642-655.
    doi: 10.1016/j.omtm.2021.04.006.
    PMID: 34095346;
    PMCID: PMC8143981.
    https://www.cell.com/molecular-therapy-family/methods/fulltext/S2329-0501(21)00072-3

    Original MATLAB source code:
    https://github.com/tamntnguyen/AAV-Triple-Transfection-Mechanistic-Model

    """

    def __init__(
        self,
        media_exchange_time: float = MEASURED_DATA.MEDIA_EXCHANGE_TIME,
        total_time: float = MEASURED_DATA.TOTAL_TIME,
        initial_concentrations=MEASURED_DATA.INITIAL_CONCENTRATIONS,
        kUptake: float = MEASURED_DATA.kUptake,
        kEscape: float = MEASURED_DATA.kEscape,
        kNuclearEntry: float = MEASURED_DATA.kNuclearEntry,
        kPlasmidDeg: float = MEASURED_DATA.kPlasmidDeg,
        kRepSyn: float = MEASURED_DATA.kRepSyn,
        kCapSyn: float = MEASURED_DATA.kCapSyn,
        kProteinDegRep: float = MEASURED_DATA.kProteinDegRep,
        kProteinDegCap: float = MEASURED_DATA.kProteinDegCap,
        kDNArep: float = MEASURED_DATA.kDNArep,
        kBindRCplasmid: float = MEASURED_DATA.kBindRCplasmid,
        kAssembly: float = MEASURED_DATA.kAssembly,
        kBindCapsid: float = MEASURED_DATA.kBindCapsid,
        kPack: float = MEASURED_DATA.kPack,
        kSecrete: float = MEASURED_DATA.kSecrete,
        ode_solver_method: str = "BDF",
    ):
        """Initialize experiment class

        Defines the class for the AAV ODE model.
        Parameters from Table 1 of the paper.
        https://www.cell.com/action/showFullTableHTML
              ?isHtml=true&tableId=tbl1&pii=S2329-0501%2821%2900072-3


        Args:
            media_exchange_time (float) : Time until media exchange
            total_time (float) : Total experiment time
            initial_concentrations (np.array) : Initial 3 plasmid concentrations
            kUptake (float) : DNA uptake from medium into cell cytosol
            kEscape (float) : endosomal escape
            kNuclearEntry (float) : nuclear entry
            kPlasmidDeg (float) : plasmid degradation
            kRepSyn (float) : Rep protein synthesis
            kCapSyn (float) : capsid protein synthesis
            kProteinDegRep (float) : Rep protein degradation
            kProteinDegCap (float) : VP degradation
            kDNArep (float) : transgene rescue and replication
            kBindRCplasmid (float) : binding of Rep protein on packaging plasmid
            kAssembly (float) : capsid assembly
            kBindCapsid (float) : binding of Rep protein on empty capsid
            kPack (float) : viral DNA encapsidation
            kSecrete (float) : capsid secretion from the nucleus to the cytosol
            ode_solver_method (str): ODE solver method (default BDF)

        """

        self.set_media_exchange_time(media_exchange_time)  # Hours
        self.set_total_time(total_time)  # Hours

        # initial plasmid input per cell for each of the 3 plasmids
        self.set_initial_concentrations(initial_concentrations)

        # Set k constants
        self.k = np.empty((len(kType),))
        self.set_k(kType.kUptake, kUptake)
        self.set_k(kType.kEscape, kEscape)
        self.set_k(kType.kNuclearEntry, kNuclearEntry)
        self.set_k(kType.kPlasmidDeg, kPlasmidDeg)
        self.set_k(kType.kRepSyn, kRepSyn)
        self.set_k(kType.kCapSyn, kCapSyn)
        self.set_k(kType.kProteinDegRep, kProteinDegRep)
        self.set_k(kType.kProteinDegCap, kProteinDegCap)
        self.set_k(kType.kDNArep, kDNArep)
        self.set_k(kType.kBindRCplasmid, kBindRCplasmid)
        self.set_k(kType.kAssembly, kAssembly)
        self.set_k(kType.kBindCapsid, kBindCapsid)
        self.set_k(kType.kPack, kPack)
        self.set_k(kType.kSecrete, kSecrete)

        # Sets print labels for outputs
        self.set_output_label()

        self.set_ode_solver_method(ode_solver_method)

    def __repr__(self) -> str:
        """Print out information about the class"""

        info = (
            click.style(
                "Ordinary Differential Equations (ODE) model "
                "of Adeno-Associated Virus (AAV) "
                "Production\n",
                fg="blue",
                bold=True,
            )
            + "Nguyen TNT, Sha S, Hong MS, Maloney AJ, Barone PW, Neufeld C,\n"
            "Wolfrum J, Springs SL, Sinskey AJ, Braatz RD.\n"
            + click.style(
                "Mechanistic model for production of "
                "recombinant adeno-associated\n"
                "virus via triple transfection of HEK293 cells.\n",
                fg="green",
                italic=True,
            )
            + "Mol Ther Methods Clin Dev. 2021 Apr 16;21:642-655.\n"
            "doi: 10.1016/j.omtm.2021.04.006.\n"
            "PMID: 34095346;\n"
            "PMCID: PMC8143981. "
            "https://www.cell.com/molecular-therapy-family/"
            "methods/fulltext/S2329-0501(21)00072-3\n"
            "\nOriginal MATLAB source code: \n"
            "https://github.com/tamntnguyen/AAV-Triple-"
            "Transfection-Mechanistic-Model\n"
        )

        return info

    #################################
    # PRIVATE FUNCTIONS
    #################################

    def _is_default_settings(self) -> bool:
        """Are the settings the default ones from the paper?

        Returns:
            bool, True if settings are the ones from the original paper
        """

        bMedia_exchange = (
            self.get_media_exchange_time() == MEASURED_DATA.MEDIA_EXCHANGE_TIME
        )
        bTotal_time = self.get_total_time() == MEASURED_DATA.TOTAL_TIME

        bInitial_concentrations = np.allclose(
            self.get_initial_concentrations(), MEASURED_DATA.INITIAL_CONCENTRATIONS
        )

        return bMedia_exchange & bTotal_time & bInitial_concentrations

    def _ode_viralProd(self, t: np.ndarray, x: np.ndarray):
        """ODE for AAV model of viral production

        Ordinary differential equations for AAV plasmid
        production.

        Nguyen TNT, Sha S, Hong MS, Maloney AJ, Barone PW, Neufeld C,
        Wolfrum J, Springs SL, Sinskey AJ, Braatz RD.
        Mechanistic model for production of recombinant adeno-associated
        virus via triple transfection of HEK293 cells.
        Mol Ther Methods Clin Dev. 2021 Apr 16;21:642-655.
        doi: 10.1016/j.omtm.2021.04.006.
        PMID: 34095346;
        PMCID: PMC8143981.
        https://www.cell.com/molecular-therapy-family/methods/
                    fulltext/S2329-0501(21)00072-3

        Original MATLAB source code:
        https://github.com/tamntnguyen/AAV-Triple-Transfection-Mechanistic-Model


        Args:
            t (np.ndarray) : Time stamps
            x (np.ndarray) : Model predictions

        Returns:
            r (list[floats]) : Model predictions (21 columns, t rows)
        """

        mu = self._calc_growth_rate(t)

        # Delivery

        dx1 = [
            -self.get_k(kType.kUptake) * x[outputTypes.pRC_extracellular],
            self.get_k(kType.kUptake) * x[outputTypes.pRC_extracellular]
            - (self.get_k(kType.kEscape) + self.get_k(kType.kPlasmidDeg) + mu)
            * x[outputTypes.pRC_endosomal],
            self.get_k(kType.kEscape) * x[outputTypes.pRC_endosomal]
            - (self.get_k(kType.kNuclearEntry) + self.get_k(kType.kPlasmidDeg) + mu)
            * x[outputTypes.pRC_cytosol],
            self.get_k(kType.kNuclearEntry) * x[outputTypes.pRC_cytosol]
            - mu * x[outputTypes.Packaging_Plasmid]
            - self.get_k(kType.kPlasmidDeg) * x[outputTypes.Packaging_Plasmid]
            - self.get_k(kType.kBindRCplasmid)
            * x[outputTypes.Packaging_Plasmid]
            * x[outputTypes.Rep_Protein],
        ]

        dx2 = [
            -self.get_k(kType.kUptake) * x[outputTypes.pVector_extracellular],
            self.get_k(kType.kUptake) * x[outputTypes.pVector_extracellular]
            - (self.get_k(kType.kEscape) + self.get_k(kType.kPlasmidDeg) + mu)
            * x[outputTypes.pVector_endosomal],
            self.get_k(kType.kEscape) * x[outputTypes.pVector_endosomal]
            - (self.get_k(kType.kNuclearEntry) + self.get_k(kType.kPlasmidDeg) + mu)
            * x[outputTypes.pVector_cytosol],
            self.get_k(kType.kNuclearEntry) * x[outputTypes.pVector_cytosol]
            - mu * x[outputTypes.Vector_Plasmid]
            - self.get_k(kType.kPlasmidDeg) * x[outputTypes.Vector_Plasmid],
        ]

        dx3 = [
            -self.get_k(kType.kUptake) * x[outputTypes.pHelper_extracellular],
            self.get_k(kType.kUptake) * x[outputTypes.pHelper_extracellular]
            - (self.get_k(kType.kEscape) + self.get_k(kType.kPlasmidDeg) + mu)
            * x[outputTypes.pHelper_endosomal],
            self.get_k(kType.kEscape) * x[outputTypes.pHelper_endosomal]
            - (self.get_k(kType.kNuclearEntry) + self.get_k(kType.kPlasmidDeg) + mu)
            * x[outputTypes.pHelper_cytosol],
            self.get_k(kType.kNuclearEntry) * x[outputTypes.pHelper_cytosol]
            - mu * x[outputTypes.Helper_Plasmid]
            - self.get_k(kType.kPlasmidDeg) * x[outputTypes.Helper_Plasmid],
        ]

        # Viral replication
        rRepSyn = (
            self.get_k(kType.kRepSyn)
            * x[outputTypes.Packaging_Plasmid]
            * x[outputTypes.Helper_Plasmid]
        )
        rCapSyn = (
            self.get_k(kType.kCapSyn)
            * x[outputTypes.Packaging_Plasmid]
            * x[outputTypes.Helper_Plasmid]
        )
        rDNArep = (
            self.get_k(kType.kDNArep)
            * x[outputTypes.Rep_Protein]
            * x[outputTypes.Vector_Plasmid]
            * x[outputTypes.Helper_Plasmid]
        )
        rBindRC = (
            self.get_k(kType.kBindRCplasmid)
            * x[outputTypes.Packaging_Plasmid]
            * x[outputTypes.Rep_Protein]
        )
        rAssembly = self.get_k(kType.kAssembly) * x[outputTypes.VP]
        rBindCapsid = (
            self.get_k(kType.kBindCapsid)
            * x[outputTypes.Rep_Protein]
            * x[outputTypes.Empty_Capsid_nucleus]
        )
        rPack = (
            self.get_k(kType.kPack)
            * x[outputTypes.Rep_bound_Empty_Plasmid]
            * x[outputTypes.vDNA]
        )

        dVP = [
            rRepSyn
            - rBindRC
            - rBindCapsid
            + rPack
            - (self.get_k(kType.kProteinDegRep) + mu) * x[outputTypes.Rep_Protein],
            rCapSyn
            - 60.0 * rAssembly
            - (self.get_k(kType.kProteinDegCap) + mu) * x[outputTypes.VP],
            rDNArep - rPack - mu * x[outputTypes.vDNA],
            rAssembly
            - (self.get_k(kType.kSecrete) + mu) * x[outputTypes.Empty_Capsid_nucleus]
            - rBindCapsid,
            rPack
            - (self.get_k(kType.kSecrete) + mu) * x[outputTypes.Full_Capsid_nucleus],
            self.get_k(kType.kSecrete) * x[outputTypes.Empty_Capsid_nucleus]
            - mu * x[outputTypes.Empty_Capsid_cytosol],
            self.get_k(kType.kSecrete) * x[outputTypes.Full_Capsid_nucleus]
            - mu * x[outputTypes.Full_Capsid_cytosol],
            rBindRC,
            rBindCapsid - rPack,
        ]

        return dx1 + dx2 + dx3 + dVP  # Vector 20 numbers (x)

    def _calc_growth_rate(self, t: np.ndarray) -> np.float32:
        """Cell growth fitted to cell growth curve as a function of time

        Args:
            t (np.ndarray): Time stamps

        Returns:
            Mean growth (mu)

        """

        # TODO: Figure out what this is
        g = [49427.9581111241, 0.0112812472592155]
        c = 1e6 + g[0] / g[1] * (1 - np.exp(-g[1] * t))
        mu = g[0] * np.exp(-g[1] * t) / c

        return mu

    def _convert_to_ml(self, t: np.ndarray, prediction: np.ndarray) -> np.ndarray:
        """Converts the predictions to mL

        Args:
            t (np.ndarray) : Time stamps
            prediction (np.ndarray) : Model predictions

        Returns:
            np.ndarray - Model predictions converted to mL
        """

        # TODO: Figure out what this is
        kCD = [49428.0668639073, 0.0112813300989750]
        CDvec = 1e6 + kCD[0] / kCD[1] * (1 - np.exp(-kCD[1] * t))

        # Convert to per mL basis
        y = np.zeros_like(prediction)

        for idx in range(np.shape(y)[0]):
            y[idx, :] = prediction[idx, :] * CDvec

        return y  # in mL

    def _fraction_of_DNA(self, prediction: np.ndarray) -> float:
        """Calculates the fraction of the DNA in the capsid

        f = sum(outputTypes.Full_Capsid_nucleus, outputTypes.Full_Capsid_cytosol) /
            sum(outputTypes.pVector_endosomal,
                outputTypes.pVector_cytosol,
                outputTypes.vDNA,
                outputTypes.Full_Capsid_nucleus,
                outputTypes.Full_Capsid_cytosol)

        Args:
            prediction (np.ndarray): Model predictions

        Returns:
            float : Percentage of viral DNA in the capsid
        """

        return (
            100.0  # Convert from fraction to percentage
            * np.sum(
                prediction[
                    [outputTypes.Full_Capsid_nucleus, outputTypes.Full_Capsid_cytosol],
                    -1,
                ]
            )
            / np.sum(
                prediction[
                    [
                        outputTypes.pVector_endosomal,
                        outputTypes.pVector_cytosol,
                        outputTypes.vDNA,
                        outputTypes.Full_Capsid_nucleus,
                        outputTypes.Full_Capsid_cytosol,
                    ],
                    -1,
                ]
            )
        )

    #################################
    # PUBLIC FUNCTIONS
    #################################

    def get_input_labels(self) -> dict:
        """Return a dictionary of all the k labels

        Returns:
            dict : A dictionary of the k constant labels
        """

        k = {
            kType.kUptake: "kUptake",
            kType.kEscape: "kEscape",
            kType.kNuclearEntry: "kNuclearEntry",
            kType.kPlasmidDeg: "kPlasmidDeg",
            kType.kRepSyn: "kRepSyn",
            kType.kCapSyn: "kCapSyn",
            kType.kProteinDegRep: "kProteinDegRep",
            kType.kProteinDegCap: "kProteinDegCap",
            kType.kDNArep: "kDNArep",
            kType.kBindRCplasmid: "kBindRCplasmid",
            kType.kAssembly: "kAssembly",
            kType.kBindCapsid: "kBindCapsid",
            kType.kPack: "kPack",
            kType.kSecrete: "kSecrete",
        }

        return k

    def get_k(self, kLabel: kType) -> float:
        """Returns one of the ODE k constants

        Args:
            kLabel(kType):  The label for the k
        Returns:
            float : The specified ODE constant.

        """
        choices = [item.value for item in kType]
        if kLabel in choices:
            return self.k[kLabel]
        else:
            raise ValueError(f"Type {kLabel} not in {choices}")

    def set_k(self, idx: kType, value: float):
        """The ODE k constants

        kUptake
        kEscape
        kNuclearEntry
        kPlasmidDeg
        kRepSyn
        kCapSyn
        kProteinDegRep [aka K_Rep_protein_degrade]
        kProteinDegCap [aka  kVP_degrade]
        kDNArep
        kBindRCplasmid
        kAssembly
        kBindCapsid
        kPack  [aka  k_DNA_pack]
        kSecrete:  Secretion constant

        Args:
            idx (kType) : k type
            value (float): value for k constant

        """

        self.k[idx] = value

    def set_output_label(
        self,
        pRC_extracellular: str = "pRC<br>extracellular",
        pRC_endosomal: str = "pRC<br>endosomal",
        pRC_cytosol: str = "pRC<br>cytosol",
        Packaging_Plasmid: str = "Packaging<br>Plasmid",  # aka pRC nucleus
        pVector_extracellular: str = "pVector<br>extracellular",
        pVector_endosomal: str = "pVector<br>endosomal",
        pVector_cytosol: str = "pVector<br>cytosol",
        Vector_Plasmid: str = "Vector<br>Plasmid",  # aka "pVector nucleus"
        pHelper_extracellular: str = "pHelper_extracellular",
        pHelper_endosomal: str = "pHelper<br>endosomal",
        pHelper_cytosol: str = "pHelper<br>cytosol",
        Helper_Plasmid: str = "Helper<br>Plasmid",  # aka "pHelper nucleus"
        Rep_Protein: str = "Rep<br>Protein",
        VP: str = "VP",  # aka "cap protein"
        vDNA: str = "vDNA",
        Empty_Capsid_nucleus: str = "Empty<br>Capsid<br>nucleus",
        Full_Capsid_nucleus: str = "Full Capsid<br>nucleus",
        Empty_Capsid_cytosol: str = "Empty Capsid<br>cytosol",
        Full_Capsid_cytosol: str = "Full Capsid<br>cytosol",
        Rep_bound_Packaging_Plasmid: str = "Rep_bound<br>Packaging<br>Plasmid",
        Rep_bound_Empty_Plasmid: str = "Rep_bound<br>Empty<br>Plasmid",
    ) -> str:
        """Set the output string

        Args:
            pRC_extracellular(str): pRC extracellular
            pRC_endosomal(str): pRC endosomal
            pRC_cytosol (str): pRC cytosol
            Packaging_Plasmid (str): Packaging Plasmid # aka pRC nucleus
            pVector_extracellular(str): pVector extracellular
            pVector_endosomal(str): pVector endosomal
            pVector_cytosol(str): pVector cytosol
            Vector_Plasmid(str): Vector Plasmid # aka "pVector nucleus"
            pHelper_extracellular(str): pHelper_extracellular
            pHelper_endosomal(str): pHelper endosomal
            pHelper_cytosol(str): pHelpercytosol
            Helper_Plasmid(str): Helper Plasmid # aka "pHelper nucleus"
            Rep_Protein(str): Rep Protein
            VP(str): VP  # aka "cap protein"
            vDNA(str): vDNA
            Empty_Capsid_nucleus(str): Empty Capsid nucleus
            Full_Capsid_nucleus(str): Full Capsid nucleus
            Empty_Capsid_cytosol(str): Empty Capsid cytosol
            Full_Capsid_cytosol(str): Full Capsid cytosol
            Rep_bound_Packaging_Plasmid(str): Rep_bound Packaging Plasmid
            Rep_bound_Empty_Plasmid(str): Rep_bound Empty Plasmid

        """

        label = {
            outputTypes.pRC_extracellular: pRC_extracellular,
            outputTypes.pRC_endosomal: pRC_endosomal,
            outputTypes.pRC_cytosol: pRC_cytosol,
            outputTypes.Packaging_Plasmid: Packaging_Plasmid,  # aka pRC nucleus
            outputTypes.pVector_extracellular: pVector_extracellular,
            outputTypes.pVector_endosomal: pVector_endosomal,
            outputTypes.pVector_cytosol: pVector_cytosol,
            outputTypes.Vector_Plasmid: Vector_Plasmid,  # aka "pVector nucleus"
            outputTypes.pHelper_extracellular: pHelper_extracellular,
            outputTypes.pHelper_endosomal: pHelper_endosomal,
            outputTypes.pHelper_cytosol: pHelper_cytosol,
            outputTypes.Helper_Plasmid: Helper_Plasmid,  # aka "pHelper nucleus"
            outputTypes.Rep_Protein: Rep_Protein,
            outputTypes.VP: VP,  # aka "cap protein"
            outputTypes.vDNA: vDNA,
            outputTypes.Empty_Capsid_nucleus: Empty_Capsid_nucleus,
            outputTypes.Full_Capsid_nucleus: Full_Capsid_nucleus,
            outputTypes.Empty_Capsid_cytosol: Empty_Capsid_cytosol,
            outputTypes.Full_Capsid_cytosol: Full_Capsid_cytosol,
            outputTypes.Rep_bound_Packaging_Plasmid: Rep_bound_Packaging_Plasmid,
            outputTypes.Rep_bound_Empty_Plasmid: Rep_bound_Empty_Plasmid,
        }

        """
        "Full cap - nucleus" + "Full cap - cytosol" = "Full Virion"

        "pVector endosomal" + "pVector cytosol" 
            + "pVector nucleus" + "vDNA" 
            + "Full cap - nucleus" + "Full cap - cytosol" 
            = "Replicated viral DNA"

        "Empty cap - nucleus" + "Full cap - nucleus" 
            + "Empty cap - cytosol" + "Full cap - cytosol" 
            = "Total capsids
        
        """

        self.output_labels = label

    def get_output_labels(self):
        """Return the string label associated with the output type

        Returns:

            list : The labels of the output predictions
        """

        labels = []
        for label in outputTypes:

            label = self.get_output_label(label).replace("<br>", " ")
            labels += [f"{label}"]

        return labels

    def get_output_label(self, index: outputTypes) -> str:
        """Return the string label associated with the output type

        Args:
            index(outputTypes): Index

        Returns:
            str: Label of this output index
        """
        choices = [item.value for item in outputTypes]
        if index in choices:
            return self.output_labels[index]
        else:
            raise ValueError(f"Type {index} not in {choices}")

    def set_media_exchange_time(
        self, media_exchange_time: float = MEASURED_DATA.MEDIA_EXCHANGE_TIME
    ):
        """Sets the time of the media exchange

        Args:
            media_exchange_time (float): The hours after transfection
            that the media is exchanged
        """
        self.MEDIA_EXCHANGE_TIME = media_exchange_time  # Hours

    def get_media_exchange_time(self) -> float:
        """Gets the time of the media exchange

        Returns:
            float : The hours after transfection that the media is exchanged
        """

        return self.MEDIA_EXCHANGE_TIME

    def set_total_time(self, total_time: float = MEASURED_DATA.TOTAL_TIME):
        """Sets the total time (hours) of the experiment

        Args:
            total_time (float): The total hours of the experiment
        """

        self.TOTAL_TIME = total_time  # Hours

    def get_total_time(self) -> float:
        """Gets the total time of the experiment

        Returns:
            float : The total time in hours of the experiment
        """

        return self.TOTAL_TIME

    def set_initial_concentrations(
        self, initial_conditions=MEASURED_DATA.INITIAL_CONCENTRATIONS
    ):
        """initial plasmid input per cell for each of the 3 plasmids"""

        self.INITIAL_PLASMID_INPUT = initial_conditions

    def get_initial_concentrations(self):
        """Get the initial concentrations of rep, cap, and helper"""

        # TODO: Is this [rep/cap plasmid, GOI plasmid, helper plasmid]?

        return self.INITIAL_PLASMID_INPUT

    def get_dna_fraction(self) -> float:
        """Get the fraction (percentage) of the DNA used by capsid

        Returns:
            float: Percentage of DNA used by capsid
        """

        if self.dna_fraction:
            return self.dna_fraction
        else:
            return None

    def get_solver_description(self) -> dict:
        """Returns a dictionary of the possible ODE solvers

        https://docs.scipy.org/doc/scipy/reference/
               generated/scipy.integrate.solve_ivp.html

        Returns:
            dict: ODE solver methods
        """
        methods = {
            "RK45": "Explicit Runge-Kutta method of order 5",
            "RK23": "Explicit Runge-Kutta method of order 3",
            "DOP853": "Explicit Runge-Kutta method of order 8",
            "Radau": "Implicit Runge-Kutta method of the Radau IIA family of order 5",
            "BDF": (
                "Implicit multi-step variable-order (1 to 5) method "
                "based on a backward differentiation formula "
                "for the derivative approximation "
            ),
            "LSODA": "Adams/BDF method with automatic stiffness detection and switching",
        }

        return methods

    def set_ode_solver_method(self, method: str = "BDF"):
        """Ordinary differential equations solver methods

        https://docs.scipy.org/doc/scipy/reference/
               generated/scipy.integrate.solve_ivp.html
        """
        self.ode_solver_method = method

    def get_ode_solver_method(self) -> str:
        """Get the ODE solver being used"""

        return self.ode_solver_method

    def run_simulation(self):
        """Initial conditions - number of plasmids

        1. Runs the ODE AAV model simulation.
        2. Starts the model based on the initial conditions.
        3. Solves the ODE from time 0 to time of media exchange.
        4. Updates the initial conditions after media exchange.
        5. Solves the ODE from time media exchange to total time.
        """

        # 2. Starts the model based on the initial conditions.
        # Run until MEDIA_EXCHANGE_TIME with initial plasmid values
        # In the paper this is 6 hours
        ic = np.zeros((len(outputTypes),))
        ic[0] = self.INITIAL_PLASMID_INPUT[0]  # pRC
        ic[4] = self.INITIAL_PLASMID_INPUT[1]  # pVector
        ic[8] = self.INITIAL_PLASMID_INPUT[2]  # pHelper

        # Time 0 to time 6 hours
        t_span = (0.0, self.get_media_exchange_time())

        # 3. Solves the ODE from time 0 to time of media exchange.
        solution = solve_ivp(
            self._ode_viralProd,
            t_span,
            ic,
            method=self.get_ode_solver_method(),
            rtol=1e-2,
            atol=1e-3,
        )

        t = solution.t
        prediction = solution.y

        if prediction is None:
            raise Exception("Parameters won't converge.")

        # Now run the remaining time after MEDIA_EXCHANGE_TIME

        # 4. Updates the initial conditions after media exchange.
        ic = prediction[:, -1]  # Just the last value of the previous prediction span

        # Media exchange
        ic[0] = 0.0  # pRC
        ic[4] = 0.0  # pVector
        ic[8] = 0.0  # pHelper

        # Time 6 hours to time 60 hours
        t_span = (self.get_media_exchange_time(), self.TOTAL_TIME)

        # 5. Solves the ODE from time media exchange to total time.
        solution = solve_ivp(
            self._ode_viralProd,
            t_span,
            ic,
            method=self.get_ode_solver_method(),
            rtol=1e-2,
            atol=1e-3,
        )

        # Append to the previous time span predictions
        t = np.append(t, solution.t)
        prediction = np.concatenate((prediction, solution.y), axis=1)

        if prediction is None:
            raise Exception("Parameters won't converge.")

        # Convert to units of mL
        prediction = self._convert_to_ml(t, prediction)

        self.dna_fraction = self._fraction_of_DNA(prediction)

        return t, prediction

    def plot_replication(self, t: np.ndarray, y: np.ndarray) -> go.Figure:
        """Plot replication curves

        This is Figure 4A from the paper.
        https://www.cell.com/molecular-therapy-family/methods/fulltext/S2329-0501(21)00072-3

        (A) Dynamics of viral GFP DNA replication.

        Args:
            t (np.ndarray):  Time stamps
            y (np.ndarray): Predictions

        Returns:
            go.Figure

        """

        idx = [
            outputTypes.pVector_endosomal,
            outputTypes.pVector_cytosol,
            outputTypes.Vector_Plasmid,
            outputTypes.vDNA,
            outputTypes.Full_Capsid_nucleus,
            outputTypes.Full_Capsid_cytosol,
        ]

        prediction = np.sum(y[idx, :], axis=0)
        max_y = np.ceil(np.max(prediction))

        fig1 = px.line(x=t, y=prediction)

        copies = MEASURED_DATA.data_Total_DNA.shape[1] - 1

        fig2 = px.scatter(
            x=np.tile(MEASURED_DATA.data_Total_DNA[:, 0], copies),
            y=np.ravel(MEASURED_DATA.data_Total_DNA[:, 1:].T),
        )

        fig1.update_traces(
            line=dict(color="black"),
            name="Replicated viral DNA - Model",
            showlegend=True,
        )

        fig3 = px.line(
            x=t,
            y=np.sum(
                y[
                    [outputTypes.Full_Capsid_nucleus, outputTypes.Full_Capsid_cytosol],
                    :,
                ],
                axis=0,
            ),
        )

        fig3.update_traces(
            line=dict(color="red"), name="Full Virion - Model", showlegend=True
        )

        if self._is_default_settings():
            fig2.update_traces(
                marker=dict(color="black", symbol="triangle-up-open", size=20),
                name="Replicated viral DNA - Measured",
                showlegend=True,
            )

            copies = MEASURED_DATA.data_full_cap.shape[1] - 1

            fig4 = px.scatter(
                x=np.tile(MEASURED_DATA.data_full_cap[:, 0], copies),
                y=np.ravel(MEASURED_DATA.data_full_cap[:, 1:].T),
            )

            fig4.update_traces(
                marker=dict(color="red", symbol="square-open", size=20),
                name="Full Virion - Measured",
                showlegend=True,
            )

            fig5 = go.Figure(data=fig1.data + fig3.data + fig2.data + fig4.data)

        else:
            fig5 = go.Figure(data=fig1.data + fig3.data)

        fig5.add_annotation(  # add a text callout with arrow
            text="Media Exchange",
            x=self.get_media_exchange_time(),
            y=max_y * 0.8,
            arrowhead=1,
            showarrow=True,
            font=dict(color="green"),
        )
        fig5.add_vline(
            x=self.get_media_exchange_time(),
            line_width=3,
            line_dash="dash",
            line_color="green",
        )

        fig5.update_layout(
            title="Viral GFP DNA replication",
            xaxis_title="Time (hours post transfection)",
            yaxis_title="AAV DNA copies/mL",
            font=dict(color="black"),
            legend=dict(orientation="h", x=0, y=-0.4, xanchor="left", yanchor="top"),
            template="gridon",
        )

        return fig5

    def plot_production(self, t: np.ndarray, y: np.ndarray):
        """Plot the predicted AAV capsid production from the model.

        This is Figure 4b from the paper.
        https://www.cell.com/molecular-therapy-family/methods/fulltext/S2329-0501(21)00072-3

        Dynamics of total and full rAAV capsid production.
        Viral production parameters were fitted to
        in-house triple-transfection experimental data.
        Measurements of total and full rAAV capsids
        are per-milliliter culture, harvested
        from both cell lysate and supernatant,
        normalized by the cell density at
        the sampling time. The dotted line
        indicates time point of media
        exchange. Initial plasmid
        input was 7.6e4/cell for each of the three plasmids.

        Args:
            t (np.ndarray):  Time stamps
            y (np.ndarray): Predictions

        Returns:
            go.Figure

        """

        """
        Don't show values < min_y
        This seems to be a bug in plotly because it should do this for us.
        """
        total_capsid = np.sum(
            y[
                [
                    outputTypes.Empty_Capsid_nucleus,
                    outputTypes.Full_Capsid_nucleus,
                    outputTypes.Empty_Capsid_cytosol,
                    outputTypes.Full_Capsid_cytosol,
                ],
                :,
            ],
            axis=0,
        )

        max_y = np.ceil(np.max(total_capsid))
        
        # No idea why min exponent doesn't work for Plotly log scale
        if max_y < 1e5:
            min_y = 1e2
        else:
            min_y = 1e7
        
        idx = total_capsid >= min_y
        total_capsid = total_capsid[idx]
        time_capsid = t[idx]

        if len(time_capsid) == 0:
            err_txt = "P1 ERROR: No virions generated<br>with these settings."
            click.secho(err_txt.replace("<br>", " "))

            fig5 = px.line(x=[0, self.get_total_time()], y=[0, 0])

            fig5.add_trace(
                go.Scatter(
                    x=[self.get_total_time() // 2],
                    y=[1e3],
                    mode="text",
                    name="P1 ERROR",
                    text=[err_txt],
                    textposition="bottom center",
                    textfont=dict(size=28),
                )
            )

        else:

            fig1 = px.line(x=time_capsid, y=total_capsid)

            fig1.update_traces(
                line=dict(color="navy"), name="Total capsid - Model", showlegend=True
            )

            total_full = np.sum(
                y[
                    [outputTypes.Full_Capsid_nucleus, outputTypes.Full_Capsid_cytosol],
                    :,
                ],
                axis=0,
            )
            idx = total_full >= min_y
            if len(idx) > 0:
                total_full = total_full[idx]
                time_capsid = t[idx]
            else:
                time_capsid = t

            if len(time_capsid) == 0:

                err_txt = "P2 ERROR: No virions generated<br>with these settings."
                click.secho(err_txt.replace("<br>", " "))
                fig5 = fig1

                fig5.add_trace(
                    go.Scatter(
                        x=[self.get_total_time() // 2],
                        y=[max(1e3, np.max(total_capsid) / 3)],
                        mode="text",
                        name="P2 ERROR",
                        text=[err_txt],
                        textposition="bottom center",
                        textfont=dict(size=28),
                    )
                )

            else:

                fig3 = px.line(x=time_capsid, y=total_full)

                fig3.update_traces(
                    line=dict(color="red"), name="Full Virion - Model", showlegend=True
                )

                if self._is_default_settings():

                    fig2 = px.scatter(
                        x=np.tile(MEASURED_DATA.data_cap[:, 0], 1),
                        y=np.ravel(MEASURED_DATA.data_cap[:, 1:].T),
                    )
                    fig2.update_traces(
                        marker=dict(color="navy", symbol="diamond-open", size=20),
                        name="Total capsid - Measured",
                        showlegend=True,
                    )

                    fig4 = px.scatter(
                        x=np.tile(MEASURED_DATA.data_full_cap[:, 0], 3),
                        y=np.ravel(MEASURED_DATA.data_full_cap[:, 1:].T),
                    )

                    fig4.update_traces(
                        marker=dict(color="red", symbol="square-open", size=20),
                        name="Full Virion - Measured",
                        showlegend=True,
                    )

                    fig5 = go.Figure(data=fig1.data + fig3.data + fig4.data + fig2.data)

                else:
                    fig5 = go.Figure(data=fig1.data + fig3.data)

        
        fig5.add_annotation(  # add a text callout with arrow
            text="Media Exchange",
            x=self.get_media_exchange_time(),
            y=max_y * 0.8,
            arrowhead=1,
            showarrow=True,
            font=dict(color="green"),
        )
        fig5.add_vline(
            x=self.get_media_exchange_time(),
            line_width=3,
            line_dash="dash",
            line_color="green",
        )

        fig5.update_layout(
            title="rAAV capsid production<br>(Total and full)",
            xaxis_title="Time (hours post transfection)",
            yaxis_title="Capsids/mL (All Compartments)",
            font=dict(color="black"),
            legend=dict(orientation="h", x=0, y=-0.4, xanchor="left", yanchor="top"),
            template="gridon",
        )
        fig5.update_yaxes(type="log")

        return fig5

    def plot_outputs(self, t: np.ndarray, y: np.ndarray):
        """
        Plot the predictions for selected model outputs.

        This is Figure 7a from the paper.
        https://www.cell.com/molecular-therapy-family/methods/fulltext/S2329-0501(21)00072-3

        Dynamic trends of the (A) species concentrations
        and (B) reaction fluxes. The values are qualitative.
        The simulation parameters were according
        to the in-house experiments described
        in Materials and methods.

        Args:
            t (np.ndarray):  Time stamps
            y (np.ndarray): Predictions

        Returns:
            go.Figure

        """

        rows = 3
        cols = 4

        plot_idx = [
            outputTypes.Packaging_Plasmid,
            outputTypes.Helper_Plasmid,
            outputTypes.Vector_Plasmid,
            outputTypes.Rep_Protein,
            outputTypes.VP,  # aka "cap protein"
            outputTypes.vDNA,
            outputTypes.Empty_Capsid_nucleus,
            outputTypes.Full_Capsid_nucleus,
            outputTypes.Empty_Capsid_cytosol,
            outputTypes.Full_Capsid_cytosol,
            outputTypes.Rep_bound_Packaging_Plasmid,
            outputTypes.Rep_bound_Empty_Plasmid,
        ]

        subplot_titles = []
        for idx in plot_idx:
            subplot_titles += [self.get_output_label(idx)]

        fig = make_subplots(
            rows=rows,
            cols=cols,
            shared_xaxes=True,
            subplot_titles=subplot_titles,
            x_title="Hours post transfection",
        )

        for row in range(rows):
            for col in range(cols):
                idx = plot_idx[col + cols * row]
                fig.add_trace(
                    go.Scatter(
                        x=t, y=y[idx], name=self.get_output_label(idx), showlegend=False
                    ),
                    row=row + 1,
                    col=col + 1,
                )
                fig.update_yaxes(visible=False)

        fig.update_layout(
            template="gridon",
            title_text="") # Blank title needed because Dash cuts off top of graph

        return fig


################################
# GLOBAL FUNCTIONS
################################


def get_version() -> str:
    """Gets the version"""

    import pkg_resources  # part of setuptools

    return pkg_resources.require(PACKAGE_NAME)[0].version


def print_version(ctx, param, value):  # pragma: no cover
    """Prints the software version"""

    if not value or ctx.resilient_parsing:
        return

    click.secho(f"ODE model of AAV production, Version {get_version()}")
    ctx.exit()


@click.command()
@click.option(
    "--version",
    "-v",
    is_flag=True,
    callback=print_version,
    expose_value=False,
    help="Display the version of this library",
    is_eager=True,
)
def main():  # pragma: no cover
    """Main function"""

    aav_model = Experiments()
    click.secho(aav_model)


if __name__ == "__main__":  # pragma: no cover

    main()
