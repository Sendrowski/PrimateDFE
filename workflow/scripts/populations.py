import hashlib

import pandas as pd


class Populations:
    """
    Population utils.
    """

    # all populations
    names = [
        'Homo_sapiens',
        'Pan_paniscus',
        'Pan_troglodytes_ellioti',
        'Pan_troglodytes_schweinfurthii',
        'Pan_troglodytes_troglodytes',
        'Pan_troglodytes_verus',
        #'Pan_troglodytes',
        #'Gorilla_gorilla',
        #'Gorilla_beringei',
        'Gorilla_beringei_beringei',
        'Gorilla_beringei_graueri',
        'Gorilla_gorilla_gorilla',
        'Gorilla_gorilla_diehli',
        'Pongo_abelii',
        'Pongo_pygmaeus',
        # 'Nomascus_leucogenys',  # only 3 individuals
        # 'Hoolock_hoolock',  # no GFF
        # 'Hylobates_lar',  # no GFF
        # 'Hylobates_moloch',  # no GFF
        # 'Macaca_fascicularis_aureus',  # no GFF
        # 'Macaca_fascicularis_fascicularis',  # no GFF
        # 'Macaca_fascicularis',  # no GFF
        'Macaca_cyclopis',
        #'Macaca_fuscata_fuscata',
        'Macaca_fuscata',
        'Macaca_mulatta',
        'Macaca_sylvanus',
        'Macaca_brunnescens',
        'Macaca_hecki',
        'Macaca_leonina',
        'Macaca_maura',
        #'Macaca_nemestrina_nemestrina',
        'Macaca_nemestrina',
        'Macaca_nigra',
        'Macaca_nigrescens',
        'Macaca_radiata',
        'Macaca_siberu',
        'Macaca_silenus',
        'Macaca_tonkeana',
        'Macaca_arctoides',
        'Macaca_assamensis',
        'Macaca_leucogenys',
        'Macaca_sinica',
        'Macaca_thibetana',
        'Papio_anubis',
        'Papio_cynocephalus_anubishybrid',
        'Papio_cynocephalus_ssp',
        #'Papio_cynocephalus', # Papio_cynocephalus_ssp + Papio_cynocephalus_anubishybrid
        'Papio_hamadryas',
        'Papio_kindae',
        # 'Papio_papio',  # folded SFS skewed
        'Papio_ursinus',
        #'Papio_ursinus_ssp',
        "Theropithecus_gelada",
        # "Chlorocebus_aethiops",  # no bcf_step1 folder for Chlorocebus
        # "Chlorocebus_cynosuros",
        # "Chlorocebus_pygerythrus",
        # "Chlorocebus_sabaeus",
        # "Chlorocebus_tantalus",
        "Erythrocebus_patas",
        "Rhinopithecus_bieti",  # folded SFS skewed but unfolded seems ok
        # "Rhinopithecus_brelichi",  # folded SFS skewed
        "Rhinopithecus_roxellana",
        # "Pygathrix_cinerea",  # folded SFS skewed
        # "Pygathrix_nemaeus",  # SFS alright, but DFE estimate very variable and genus appears problematic
        # "Pygathrix_nigripes",  # folded SFS skewed
        # "Semnopithecus_entellus",  # unfolded SFS skewed, also folded SFS to some extent
        "Semnopithecus_hypoleucos",
        "Trachypithecus_francoisi",
        # "Trachypithecus_geei",  # folded SFS skewed
        # "Trachypithecus_phayrei",  # folded SFS skewed
        "Trachypithecus_poliocephalus",  # folded SFS skewed but unfolded seems ok
        "Callithrix_jacchus",
        "Saguinus_oedipus",
        "Saguinus_midas",
        "Leontocebus_fuscicollis",
        "Aotus_azarai",
        # "Aotus_griseimembra",  # folded SFS skewed
        "Aotus_nancymaae",
        # "Aotus_vociferans",  # folded SFS skewed
        "Cebus_albifrons",
        "Cebus_unicolor",
        "Sapajus_apella_macrocephalus",
        # "Saimiri_boliviensis_boliviensis",  # folded SFS skewed
        "Saimiri_cassiquiarensis",
        "Saimiri_ustus",
        # "Alouatta_macconnelli", # reference has no GFF
        # "Alouatta_palliata",
        # "Alouatta_seniculus_puruensis",
        # "Alouatta_seniculus_ssp",
        # "Ateles_chamek",
        "Plecturocebus_cupreus",
        "Cheracebus_lucifer",
        "Cheracebus_lugens",
        # "Cacajao_calvus_rubicundus", # reference has no GFF
        # "Cacajao_calvus_ssp",
        # "Chiropotes_albinasus",
        # "Chiropotes_israelita",
        # "Chiropotes_sagulatus",
        # "Pithecia_pithecia"
    ]

    groups = {}

    groups["great_apes"] = [
        'Homo_sapiens',
        #'Pan_troglodytes',
        'Pan_paniscus',
        'Pan_troglodytes_ellioti',
        'Pan_troglodytes_schweinfurthii',
        'Pan_troglodytes_troglodytes',
        'Pan_troglodytes_verus',
        #'Gorilla_gorilla',
        #'Gorilla_beringei',
        'Gorilla_beringei_beringei',
        'Gorilla_beringei_graueri',
        'Gorilla_gorilla_gorilla',
        'Gorilla_gorilla_diehli',
        'Pongo_abelii',
        'Pongo_pygmaeus'
    ]

    groups["Castellano_data"] = [
        "Homo_sapiens",
        "Pan_paniscus",
        "Pan_troglodytes_troglodytes",
        "Pan_troglodytes_schweinfurthii",
        "Pan_troglodytes_verus",
        "Pan_troglodytes_ellioti",
        "Pongo_pygmaeus",
        "Pongo_abelii",
        "Gorilla_gorilla_gorilla",
    ]

    groups["Castellano"] = [
        "Homo_sapiens",
        "Pan_paniscus",
        "Pan_troglodytes_troglodytes",
        "Pan_troglodytes_schweinfurthii",
        "Pan_troglodytes_verus",
        "Pan_troglodytes_ellioti",
        "Pongo_pygmaeus",
        "Pongo_abelii",
        "Gorilla_gorilla_gorilla",
    ]

    groups["macaques"] = [p for p in names if p.startswith('Macaca_')]

    groups["baboons"] = [p for p in names if p.startswith('Papio_')]

    groups["cheek_pouch_monkeys"] = (
            groups["macaques"] +
            groups["baboons"] +
            # [p for p in populations if p.startswith('Chlorocebus_')] +
            ["Theropithecus_gelada", "Erythrocebus_patas"]
    )

    groups["leaf_eating_monkeys"] = [
        p for p in names
        if (
                p.startswith("Rhinopithecus_")
                or p.startswith("Pygathrix_")
                or p.startswith("Semnopithecus_")
                or p.startswith("Trachypithecus_")
        )
    ]

    groups["old_world_monkeys"] = groups["cheek_pouch_monkeys"] + groups["leaf_eating_monkeys"]

    groups['new_world_monkeys'] = [
        "Callithrix_jacchus",
        "Saguinus_midas",
        "Leontocebus_fuscicollis",
        "Aotus_azarai",
        # "Aotus_griseimembra",
        "Aotus_nancymaae",
        # "Aotus_vociferans",
        "Cebus_albifrons",
        "Cebus_unicolor",
        "Sapajus_apella_macrocephalus",
        # "Saimiri_boliviensis_boliviensis",
        "Saimiri_cassiquiarensis",
        "Saimiri_ustus",
        # "Alouatta_macconnelli",
        # "Alouatta_palliata",
        # "Alouatta_seniculus_puruensis",
        # "Alouatta_seniculus_ssp",
        # "Ateles_chamek",
        "Plecturocebus_cupreus",
        "Cheracebus_lucifer",
        "Cheracebus_lugens",
        # "Cacajao_calvus_rubicundus",
        # "Cacajao_calvus_ssp",
        # "Chiropotes_albinasus",
        # "Chiropotes_israelita",
        # "Chiropotes_sagulatus",
        # "Pithecia_pithecia"
    ]

    groups["catarrhini"] = groups["great_apes"] + groups["old_world_monkeys"]

    groups["all"] = names

    names_castellano = [
        "Homo_sapiens",  # human
        "Pan_paniscus",  # bonobo
        "Pan_troglodytes_troglodytes",  # central_chimp
        "Pan_troglodytes_schweinfurthii",  # eastern_chimp
        "Pan_troglodytes_verus",  # western_chimp
        "Pan_troglodytes_ellioti",  # NC_chimp
        "Pongo_pygmaeus",  # bornean_orang
        "Pongo_abelii",  # sumatran_orang
        "Gorilla_gorilla_gorilla",  # western_lowland_gorilla
    ]

    names_castellano_english = [
        'human',
        'bonobo',
        'central_chimp',
        'eastern_chimp',
        'western_chimp',
        'NC_chimp',
        'bornean_orang',
        'sumatran_orang',
        'western_lowland_gorilla',
    ]

    # outgroups for each ingroup species
    ingroup_outgroups = dict(
        Homo_sapiens=['Pan_troglodytes', 'Gorilla_gorilla', 'Pongo_abelii'],
        Pan_troglodytes=['Homo_sapiens', 'Gorilla_gorilla', 'Pongo_abelii'],
        Pan_paniscus=['Homo_sapiens', 'Gorilla_gorilla', 'Pongo_abelii'],
        Gorilla_gorilla=['Homo_sapiens', 'Pongo_abelii'],
        Gorilla_beringei=['Homo_sapiens', 'Pongo_abelii'],
        Pongo_abelii=['Homo_sapiens', 'Nomascus_leucogenys'],
        Pongo_pygmaeus=['Homo_sapiens', 'Nomascus_leucogenys'],
        Nomascus_leucogenys=['Hylobates_pileatus', 'Homo_sapiens'],
        Hoolock_hoolock=['Hylobates_pileatus', 'Nomascus_leucogenys', 'Homo_sapiens'],
        Hoolock_leuconedys=['Hylobates_pileatus', 'Nomascus_leucogenys', 'Homo_sapiens'],
        Hylobates_lar=['Hylobates_pileatus', 'Nomascus_leucogenys', 'Homo_sapiens'],
        Hylobates_moloch=['Hylobates_pileatus', 'Nomascus_leucogenys', 'Homo_sapiens'],
        Macaca_cyclopis=['Macaca_thibetana', 'Macaca_nemestrina', 'Papio_anubis'],
        Macaca_fuscata_fuscata=['Macaca_thibetana', 'Macaca_nemestrina', 'Papio_anubis'],
        Macaca_fuscata=['Macaca_thibetana', 'Macaca_nemestrina', 'Papio_anubis'],
        Macaca_mulatta=['Macaca_thibetana', 'Macaca_nemestrina', 'Papio_anubis'],
        Macaca_sylvanus=['Macaca_thibetana', 'Macaca_nemestrina', 'Papio_anubis'],
        Macaca_hecki=['Macaca_nemestrina', 'Macaca_mulatta', 'Papio_anubis'],
        Macaca_leonina=['Macaca_nemestrina', 'Macaca_mulatta', 'Papio_anubis'],
        Macaca_maura=['Macaca_nemestrina', 'Macaca_mulatta', 'Papio_anubis'],
        Macaca_nemestrina_nemestrina=['Macaca_leonina', 'Macaca_mulatta', 'Papio_anubis'],
        Macaca_nemestrina=['Macaca_leonina', 'Macaca_mulatta', 'Papio_anubis'],
        Macaca_brunnescens=['Macaca_leonina', 'Macaca_mulatta', 'Papio_anubis'],
        Macaca_nigrescens=['Macaca_leonina', 'Macaca_mulatta', 'Papio_anubis'],
        Macaca_siberu=['Macaca_leonina', 'Macaca_mulatta', 'Papio_anubis'],
        Macaca_nigra=['Macaca_nemestrina', 'Macaca_mulatta', 'Papio_anubis'],
        Macaca_radiata=['Macaca_mulatta', 'Macaca_nemestrina', 'Papio_anubis'],
        Macaca_thibetana=['Macaca_mulatta', 'Macaca_nemestrina', 'Papio_anubis'],
        Macaca_sinica=['Macaca_mulatta', 'Macaca_nemestrina', 'Papio_anubis'],
        Macaca_silenus=['Macaca_nemestrina', 'Macaca_mulatta', 'Papio_anubis'],
        Macaca_tonkeana=['Macaca_nemestrina', 'Macaca_mulatta', 'Papio_anubis'],
        Macaca_arctoides=['Macaca_mulatta', 'Macaca_nemestrina', 'Papio_anubis'],
        Macaca_assamensis=['Macaca_mulatta', 'Macaca_nemestrina', 'Papio_anubis'],
        Macaca_leucogenys=['Macaca_mulatta', 'Macaca_nemestrina', 'Papio_anubis'],
        Papio_anubis=['Papio_hamadryas', 'Papio_kindae', 'Macaca_mulatta'],
        Papio_cynocephalus=['Papio_ursinus', 'Papio_anubis', 'Macaca_mulatta'],
        Papio_hamadryas=['Papio_anubis', 'Papio_kindae', 'Macaca_mulatta'],
        Papio_kindae=['Papio_cynocephalus', 'Papio_anubis', 'Macaca_mulatta'],
        Papio_papio=['Papio_hamadryas', 'Papio_kindae', 'Macaca_mulatta'],
        Papio_ursinus=['Papio_cynocephalus', 'Papio_anubis', 'Macaca_mulatta'],
        Theropithecus_gelada=['Papio_anubis', 'Macaca_mulatta', 'Allenopithecus_nigroviridis'],
        Chlorocebus_aethiops=['Chlorocebus_pygerythrus', 'Allenopithecus_nigroviridis', 'Macaca_mulatta'],
        Chlorocebus_cynosuros=['Chlorocebus_sabaeus', 'Allenopithecus_nigroviridis', 'Macaca_mulatta'],
        Chlorocebus_pygerythrus=['Chlorocebus_sabaeus', 'Allenopithecus_nigroviridis', 'Macaca_mulatta'],
        Chlorocebus_sabaeus=['Chlorocebus_tantalus', 'Allenopithecus_nigroviridis', 'Macaca_mulatta'],
        Chlorocebus_tantalus=['Chlorocebus_sabaeus', 'Allenopithecus_nigroviridis', 'Macaca_mulatta'],
        Erythrocebus_patas=['Chlorocebus_sabaeus', 'Papio_anubis', 'Macaca_mulatta'],
        Rhinopithecus_bieti=['Rhinopithecus_roxellana', 'Trachypithecus_francoisi', 'Colobus_guereza'],
        Rhinopithecus_brelichi=['Rhinopithecus_roxellana', 'Trachypithecus_francoisi', 'Colobus_guereza'],
        Rhinopithecus_roxellana=['Pygathrix_cinerea', 'Trachypithecus_francoisi', 'Colobus_guereza'],
        Pygathrix_cinerea=['Rhinopithecus_roxellana', 'Trachypithecus_francoisi', 'Colobus_guereza'],
        Pygathrix_nemaeus=['Rhinopithecus_bieti', 'Trachypithecus_francoisi', 'Colobus_guereza'],
        Pygathrix_nigripes=['Rhinopithecus_brelichi', 'Trachypithecus_francoisi', 'Colobus_guereza'],
        Semnopithecus_entellus=['Trachypithecus_francoisi', 'Rhinopithecus_roxellana', 'Colobus_guereza'],
        Semnopithecus_hypoleucos=['Trachypithecus_francoisi', 'Rhinopithecus_roxellana', 'Colobus_guereza'],
        Trachypithecus_francoisi=['Semnopithecus_entellus', 'Trachypithecus_francoisi', 'Colobus_guereza'],
        Trachypithecus_geei=['Trachypithecus_phayrei', 'Semnopithecus_entellus', 'Colobus_guereza'],
        Trachypithecus_phayrei=['Trachypithecus_geei', 'Semnopithecus_entellus', 'Colobus_guereza'],
        Trachypithecus_poliocephalus=['Trachypithecus_geei', 'Semnopithecus_entellus', 'Colobus_guereza'],
        Callithrix_jacchus=['Saguinus_oedipus', 'Sapajus_apella'],
        Saguinus_oedipus=['Saguinus_midas', 'Callithrix_jacchus'],
        Saguinus_midas=['Saguinus_oedipus', 'Callithrix_jacchus'],
        Leontocebus_fuscicollis=['Saguinus_oedipus', 'Callithrix_jacchus', 'Aotus_nancymaae'],
        Aotus_azarai=['Saguinus_oedipus', 'Sapajus_apella'],
        Aotus_griseimembra=['Saguinus_oedipus', 'Sapajus_apella'],
        Aotus_nancymaae=['Saguinus_oedipus', 'Sapajus_apella'],
        Aotus_vociferans=['Saguinus_oedipus', 'Sapajus_apella'],
        Cebus_albifrons=['Sapajus_apella', 'Saimiri_boliviensis'],
        Cebus_unicolor=['Sapajus_apella', 'Saimiri_boliviensis'],
        Sapajus_apella=['Saimiri_boliviensis', 'Saguinus_oedipus'],
        Saimiri_boliviensis=['Sapajus_apella', 'Saguinus_oedipus'],
        Saimiri_cassiquiarensis=['Sapajus_apella', 'Saguinus_oedipus'],
        Saimiri_ustus=['Sapajus_apella', 'Saguinus_oedipus'],
        Alouatta_macconnelli=['Ateles_hybridus', 'Sapajus_apella'],
        Alouatta_palliata=['Ateles_hybridus', 'Sapajus_apella'],
        Alouatta_seniculus=['Ateles_hybridus', 'Sapajus_apella'],
        Ateles_chamek=['Ateles_hybridus', 'Sapajus_apella'],
        Plecturocebus_cupreus=['Pithecia_pithecia', 'Sapajus_apella'],
        Cheracebus_lucifer=['Pithecia_pithecia', 'Sapajus_apella'],
        Cheracebus_lugens=['Pithecia_pithecia', 'Sapajus_apella'],
        Cacajao_calvus=['Pithecia_pithecia', 'Plecturocebus_cupreus'],
        Chiropotes_albinasus=['Pithecia_pithecia', 'Plecturocebus_cupreus'],
        Chiropotes_israelita=['Pithecia_pithecia', 'Plecturocebus_cupreus'],
        Chiropotes_sagulatus=['Pithecia_pithecia', 'Plecturocebus_cupreus'],
        Pithecia_pithecia=['Plecturocebus_cupreus', 'Sapajus_apella'],
    )

    # human chromosomes
    human_chrs = [f"chr{chr}" for chr in range(1, 23)] + ['chrX']

    dfe_vs_ne_stats = ['S_d', 's_d', 'range_S_inf_-10', 'range_S_-10_-1', 'range_S_-1_0']

    @staticmethod
    def get_group_from_pop(pop: str) -> str:
        """
        Get the group name from a population name.
        """
        for group, pops in Populations.groups.items():
            if any(pop in p for p in pops):
                return group

        return "unknown"

    @staticmethod
    def get_pops(n: int, group: str = None, path="results/stats/sample_names/counts/population_counts.csv") -> list:
        """
        Get populations with at least n haplotypes.
        """
        counts = pd.read_csv(path)

        names = counts[counts.n_individuals * 2 >= int(n)].population.tolist()

        if group:
            names = [n for n in names if n in Populations.groups[group]]

        return names

    @staticmethod
    def get_individuals_from_pop(pop: str) -> list:
        """
        Get individual ids from population name.

        :param pop: Population name or species name.
        """
        genus = pop.split('_')[0]

        df = pd.read_csv(f"resources/metadata/{genus}_individuals.csv", sep="\t")

        return df[df["BAM_FOLDER"].str.contains(pop)].GVCF_ID.tolist()

    @staticmethod
    def get_label_rank(label: str) -> int:
        """
        Get the order of labels as they first appear in the list.
        """
        rank = dict(
            humans=0,
            great_apes=1,
            macaques=2,
            baboons=3,
            cheek_pouch_monkeys=4,
            leaf_eating_monkeys=5
        )

        if label in rank:
            return rank[label]

        return int(hashlib.md5(label.encode()).hexdigest(), 16) % 100

    @staticmethod
    def get_label_color_map(labels: list) -> dict:
        """
        Return a dict mapping label -> color.
        Soft, manually chosen colors (not too strong),
        stable across all figures.
        """
        ordered_labels = sorted(
            set(labels),
            key=Populations.get_label_rank,
        )

        return {
            label: Populations.get_color(label) for label in ordered_labels
        }

    @staticmethod
    def get_color(label: str) -> str:
        """
        Get color for a given population label.
        """
        palette = {
            "great_apes": "#4C72B0",  # muted blue
            "leaf_eating_monkeys": "#55A868",  # green
            "baboons": "#8172B2",  # purple
            "macaques": "#CCB974",  # ochre
            "cheek_pouch_monkeys": "#DD8452",  # muted orange
        }

        if label in palette:
            return palette[label]

        group = Populations.get_group_from_pop(label)
        if group in palette:
            return palette[group]

        return "#999999"

    @staticmethod
    def label_to_text(label: str) -> str:
        """
        Convert population label to more readable text.
        """
        return label.replace('_', ' ')

    @staticmethod
    def to_species(name: str) -> str:
        """
        Remove subspecies or population suffix if present.
        """
        return name.rsplit("_", 1)[0] if name.count("_") >= 2 else name
