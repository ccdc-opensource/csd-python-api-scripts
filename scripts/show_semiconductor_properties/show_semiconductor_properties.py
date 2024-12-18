#
# This script can be used for any purpose without limitation subject to the
# conditions at http://www.ccdc.cam.ac.uk/Community/Pages/Licences/v2.aspx
#
# This permission notice and the following statement of attribution must be
# included in all copies or substantial portions of this script.
#
# 2023-08-18 Created by Chris Kingsbury, the Cambridge Crystallographic Data Centre
# ORCID 0000-0002-4694-5566
#
#
from ccdc.utilities import ApplicationInterface
import ccdc.io
from jinja2 import Template
from pathlib import Path
import os
from io import BytesIO
import base64
import json
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

README_LINK = f"""{Path(os.path.dirname(__file__)) / "ReadMe.md"}"""
README_LINK = "https://downloads.ccdc.cam.ac.uk/documentation/API/descriptive_docs/predicted_properties.html"

default_settings = {}


def plot_hist(descs_data, astype="fig"):
    data_hlt = {
        "singlet_state_1_energy": descs_data.singlet_state_1_energy,
        "singlet_state_1_oscillator_strength": descs_data.singlet_state_1_oscillator_strength,
        "singlet_state_2_energy": descs_data.singlet_state_2_energy,
        "singlet_state_2_oscillator_strength": descs_data.singlet_state_2_oscillator_strength,
        "triplet_state_1_energy": descs_data.triplet_state_1_energy,
        "triplet_state_2_energy": descs_data.triplet_state_2_energy,
        "homo_lumo_gap": descs_data.homo_lumo_gap,
        "transfer_integral": descs_data.transfer_integral,
        "hole_reorganization_energy": descs_data.hole_reorganization_energy,
        "dynamic_disorder": descs_data.dynamic_disorder,
    }
#    output1 = open("output1.txt", "a")
#    output1.write("%s" % (str(data_hlt{"singlet_state_1_energy"})))
    hist_data = json.load(open(Path(__file__).parent / "hist_data.json", "r"))
    fig, axs = plt.subplots(nrows=2, ncols=5, figsize=(15, 15))

    for ix, (key, data_red) in enumerate(hist_data.items()):
        row, col = int(np.floor(ix / 5)), ix % 5
        ax = axs[row][col]
        xs = data_red["x"]
        ys = data_red["y"]

        ax.stairs(xs, ys, orientation="horizontal", fill=True)
        ax.title.set_text(data_red["axis_label"])

        arrow_place = data_hlt.get(data_red["name"], np.nan)
        if arrow_place is None:
            arrow_place = np.nan
        try:
            arrow_anchor = (0, [x for x in ys if (x < arrow_place)][-1])
            height = ys[1] - ys[0]
            width = [x for x, y in zip(xs, ys) if (y < arrow_place)][-1]
            axs[row][col].add_patch(
                Rectangle(
                    xy=arrow_anchor,
                    width=width,
                    height=height,
                    facecolor="red",
                    fill=True,
                )
            )
        except IndexError:
            pass

    if astype == "buf":
        buf = BytesIO()
        fig.savefig(buf, format="png", dpi=600, backend="Agg")
        image_base64 = (
            base64.b64encode(buf.getvalue()).decode("utf-8").replace("\n", "")
        )
        buf.close()
        return image_base64

    else:
        return fig


def write_descs_report(settings=default_settings):
    interface = ApplicationInterface(parse_commandline=False)
    interface.parse_commandline()
    semiconductor_entry_reader = ccdc.io.EntryReader("CSD")
    entry = semiconductor_entry_reader.entry(interface.identifier)
    entry = interface.current_entry
    if (entry.predicted_properties is None):
        interface.write_report(title="Data not found", content="No Predicted Property Data Found For " + entry.identifier)
        return None
    else:
        properties = entry.predicted_properties
    if (properties.semiconductor_properties is None):
        interface.write_report(title="Data not found", content="No Semiconductor Data Found For " + entry.identifier)
        return None

    descs_data = properties.semiconductor_properties
    with open(interface.output_html_file, "w") as report:

        tl = Template(
            open(
                Path(__file__).parent / "semiconductor_template.html",
                "r",
            ).read()
        )
        report.write(
            tl.render(
                ident=interface.identifier, data=descs_data, readme_link=README_LINK, image=plot_hist(descs_data, "buf")
            )
        )


if __name__ == "__main__":
    write_descs_report()
