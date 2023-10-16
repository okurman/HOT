from os.path import join
import pandas as pd

from overbinders.data_prep.basic import load_metadata

data_path = "/panfs/pan1/devdcode/sanjar/overbinders"


def extract_LLPSDB_proteins():

    a = load_metadata()

    llps_file = join(data_path, "LLPS_proteins/LLPSDB/LLPS.xls")
    llps = pd.read_excel(llps_file)

    llps = llps[llps["Phase separation"] == 'Yes']
    llps_protein_ids = set(llps["Protein ID"].values)

    p_df = pd.read_excel(join(data_path, "LLPS_proteins/LLPSDB/protein.xls"))
    p_df = p_df.loc[(p_df["PID"].isin(llps_protein_ids)) & (p_df["Species"] == "Homo sapiens")]

    pnas = p_df["Protein name abbreviation"].values

    all_tfs = list(a["HepG2"].keys()) + list(a["K562"].keys())

    for pna in pnas:

        matched_tfs = [tf for tf in all_tfs if pna.lower() in tf.lower()]
        if matched_tfs:
            print(pna, matched_tfs)


def get_LLSPDB_proteins():

    # the list is selected using the function above
    tf_list = ['FUS', 'PTBP1', 'HNRNPH1', 'BRD4', 'MED1', 'MED13', 'EWSR1', 'TAF15', 'PRPF4', 'MYC', 'RARA', 'GATA2',
               'ERF', 'CBX2']

    return tf_list


def get_mierlo_etal_proteins():

    # list is obtained from https://www.sciencedirect.com/science/article/pii/S2211124721000188#mmc3

    tfs = ["AKAP95", "ATX2", "BRD3", "BRD4", "C9ORF72", "CBX2", "CBX5", "CDCA8", "COIL", "CPEB2", "CPSF6", "DCP1A",
           "DDX3X", "DDX4", "DYRK3", "EIF4G1", "EIF4H", "ELN", "EWSR1", "FBL", "FMR1", "FUS", "G3BP1", "G3BP2", "GLI3",
           "GM130", "GRB2", "HNRNPA0", "HNRNPA1", "HNRNPA1L2", "HNRNPA2B1", "HNRNPA3", "HNRNPDL", "HNRNPH1", "HNRNPH2",
           "HNRNPH3", "IAPP", "IGF2BP2", "LEMD2", "MAPT", "MATR3", "MECP2", "MED1", "NCL", "NELFE", "NFE2L2", "NONO",
           "NPM1", "NR3C1", "NUP98", "OCT4", "PABD1", "PLK4", "PML", "PRKAR1A", "PSPC1", "RAD23B", "RAD52", "RBM14",
           "RBM3 ", "RXRG", "SAFB", "SFPQ", "SMN1 ", "SNRNP70", "SOD1", "SOS1", "SP100", "SQSTM1", "SRSF1", "SRSF2",
           "STAU", "SURF6", "TACC3", "TAF15", "TARDBP", "TAZ", "TERF2", "TIA1", "TIAL1", "TPX2", "TRNP1", "UBQLN2",
           "YAP1", "YTHDF1", "YTHDF2", "YTHDF3"]

    return tfs


def get_llps_tfs():

    tfs = set(get_LLSPDB_proteins() + get_mierlo_etal_proteins())

    return tfs
