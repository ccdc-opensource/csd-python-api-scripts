#
# This script can be used for any purpose without limitation subject to the
# conditions at http://www.ccdc.cam.ac.uk/Community/Pages/Licences/v2.aspx
#
# This permission notice and the following statement of attribution must be
# included in all copies or substantial portions of this script.
#
# The following line states a licence feature that is required to show this script in Mercury and Hermes script menus.
# Created 18/08/2024 by Alex Moldovan (https://orcid.org/0000-0003-2776-3879)


import os
import sys
import tkinter as tk
from tkinter import ttk, messagebox, filedialog

from ccdc.utilities import ApplicationInterface

from _surface_charge_calculator import SurfaceChargeController


class SurfaceChargeGUI:
    def __init__(self, initial_file_path=None):
        self.root = tk.Tk()
        self.root.title("Surface Charge Calculator")
        try:
            photo = tk.PhotoImage(file=os.path.join(os.path.dirname(__file__), 'assets/csd-python-api-logo.png'))
            self.root.wm_iconphoto(False, photo)
        except FileNotFoundError:
            print("Could not find icon file for app.")
        except Exception as e:
            print("Unable to load icon")
            print(e)  # This doesn't seem to work with X11 port forwarding 🤷‍♀️
        # Disable window resizing
        self.root.resizable(False, False)

        self.initial_file_path = initial_file_path
        self.create_string_file_inputs()
        self.create_input_fields()
        self.create_buttons()
        self.create_treeview()
        self.create_directory_selection()
        self.configure_grid()  # Ensure grid configuration
        if self.initial_file_path:
            self.handle_initial_file_path(self.initial_file_path)

    def handle_initial_file_path(self, file_path):
        """Handles the initial file path by disabling the input fields and setting the file path."""
        self.file_var.set(file_path)  # Set the provided file path
        self.string_var.set("")  # Clear the string input

        # Disable the input fields
        self.string_entry.config(state='disabled')
        self.file_entry.config(state='readonly')
        self.browse_button.config(state='disabled')

    def configure_grid(self):
        self.root.grid_rowconfigure(8, weight=1)
        self.root.grid_rowconfigure(9, weight=0)
        self.root.grid_rowconfigure(10, weight=0)

        self.root.grid_columnconfigure(0, weight=1)
        self.root.grid_columnconfigure(1, weight=1)
        self.root.grid_columnconfigure(2, weight=1)
        self.root.grid_columnconfigure(3, weight=1)
        self.root.grid_columnconfigure(4, weight=1)
        self.root.grid_columnconfigure(5, weight=1)
        self.root.grid_columnconfigure(6, weight=1)
        self.root.grid_columnconfigure(7, weight=1)

    def create_string_file_inputs(self):
        tk.Label(self.root, text="Structure").grid(row=0, column=0, columnspan=2, sticky='w')

        tk.Label(self.root, text="Refcode:").grid(row=1, column=0, padx=5, pady=5, sticky='e')
        self.string_var = tk.StringVar()
        self.string_entry = tk.Entry(self.root, textvariable=self.string_var, validate="key",
                                     validatecommand=(self.root.register(self.on_string_input), "%P"))
        self.string_entry.grid(row=1, column=1, padx=5, pady=5, columnspan=2, sticky='ew')

        tk.Label(self.root, text="Select File:").grid(row=2, column=0, padx=5, pady=5, sticky='e')
        self.file_var = tk.StringVar()
        self.file_entry = tk.Entry(self.root, textvariable=self.file_var, state='readonly')
        self.file_entry.grid(row=2, column=1, padx=5, pady=5, columnspan=2, sticky='ew')
        self.browse_button = tk.Button(self.root, text="Browse", command=self.browse_file)
        self.browse_button.grid(row=2, column=3, padx=5, pady=5, sticky='ew')

    def on_string_input(self, input_value):
        if input_value.strip():
            self.browse_button.config(state='disabled')
        else:
            self.browse_button.config(state='normal')
        return True

    def create_input_fields(self):
        tk.Label(self.root, text="Select hkl and offset").grid(row=3, column=0, columnspan=2, sticky='w')

        input_frame = tk.Frame(self.root)
        input_frame.grid(row=4, column=0, columnspan=8, padx=5, pady=5, sticky='ew')

        input_frame.grid_columnconfigure(0, weight=1)
        input_frame.grid_columnconfigure(1, weight=1)
        input_frame.grid_columnconfigure(2, weight=1)
        input_frame.grid_columnconfigure(3, weight=1)
        input_frame.grid_columnconfigure(4, weight=1)
        input_frame.grid_columnconfigure(5, weight=1)
        input_frame.grid_columnconfigure(6, weight=1)
        input_frame.grid_columnconfigure(7, weight=1)

        tk.Label(input_frame, text="h:").grid(row=0, column=0, padx=2, pady=5, sticky='e')
        tk.Label(input_frame, text="k:").grid(row=0, column=2, padx=2, pady=5, sticky='e')
        tk.Label(input_frame, text="l:").grid(row=0, column=4, padx=2, pady=5, sticky='e')
        tk.Label(input_frame, text="offset:").grid(row=0, column=6, padx=2, pady=5, sticky='e')

        self.h_var = tk.IntVar()
        self.spin_h = tk.Spinbox(input_frame, from_=-9, to=9, width=2, textvariable=self.h_var)
        self.spin_h.grid(row=0, column=1, padx=2, pady=5, sticky='ew')

        self.k_var = tk.IntVar()
        self.spin_k = tk.Spinbox(input_frame, from_=-9, to=9, width=2, textvariable=self.k_var)
        self.spin_k.grid(row=0, column=3, padx=2, pady=5, sticky='ew')

        self.l_var = tk.IntVar()
        self.spin_z = tk.Spinbox(input_frame, from_=-9, to=9, width=2, textvariable=self.l_var)
        self.spin_z.grid(row=0, column=5, padx=2, pady=5, sticky='ew')

        self.offset_var = tk.DoubleVar()
        self.entry_offset = tk.Entry(input_frame, width=10, textvariable=self.offset_var)
        self.entry_offset.grid(row=0, column=7, padx=2, pady=5, sticky='ew')

    def create_buttons(self):
        self.add_button = tk.Button(self.root, text="Add Surface", command=self.add_combination)
        self.add_button.grid(row=5, column=0, columnspan=2, pady=10, sticky='ew')

        self.delete_button = tk.Button(self.root, text="Delete Selected", command=self.delete_combination)
        self.delete_button.grid(row=5, column=2, pady=5, sticky='ew')

        self.reset_button = tk.Button(self.root, text="Reset Fields", command=self.reset_fields)
        self.reset_button.grid(row=5, column=3, pady=5, sticky='ew')

        self.create_directory_selection()

    def create_directory_selection(self):
        tk.Label(self.root, text="Output Directory:").grid(row=9, column=0, padx=5, pady=5, sticky='e')

        self.dir_var = tk.StringVar(value=os.getcwd())  # Default to current working directory
        self.dir_entry = tk.Entry(self.root, textvariable=self.dir_var, state='readonly', width=50)
        self.dir_entry.grid(row=9, column=1, padx=5, pady=5, columnspan=3, sticky='ew')

        self.browse_dir_button = tk.Button(self.root, text="Browse", command=self.select_directory)
        self.browse_dir_button.grid(row=9, column=4, padx=5, pady=5, sticky='ew')

        self.calculate_button = tk.Button(self.root, text="Calculate", command=self.calculate)
        self.calculate_button.grid(row=10, column=0, columnspan=5, pady=10, sticky='ew')

    def select_directory(self):
        selected_dir = filedialog.askdirectory(initialdir=self.dir_var.get(), title="Select Output Directory")
        if selected_dir:
            self.dir_var.set(selected_dir)

    def create_treeview(self):

        tk.Label(self.root, text="Current Selections").grid(row=7, column=0, padx=5, pady=5, columnspan=8,
                                                            sticky='w')
        self.combination_tree = ttk.Treeview(self.root, columns=("h", "k", "l", "Offset"), show='headings')
        self.combination_tree.grid(row=8, column=0, columnspan=8, padx=10, pady=10, sticky='nsew')

        self.combination_tree.heading("h", text="h")
        self.combination_tree.heading("k", text="k")
        self.combination_tree.heading("l", text="l")
        self.combination_tree.heading("Offset", text="Offset")

        self.combination_tree.column("h", width=50, anchor=tk.CENTER)
        self.combination_tree.column("k", width=50, anchor=tk.CENTER)
        self.combination_tree.column("l", width=50, anchor=tk.CENTER)
        self.combination_tree.column("Offset", width=100, anchor=tk.CENTER)

    def browse_file(self):
        file_path = filedialog.askopenfilename(filetypes=[("mol2 files", "*.mol2")])
        if file_path:
            self.file_var.set(file_path)

    def add_combination(self):
        try:
            h = self.h_var.get()
            k = self.k_var.get()
            l = self.l_var.get()
            if (h, k, l) == (0, 0, 0):
                messagebox.showerror("Invalid input", "Please enter valid integers for h, k, l and a float for offset.")
                return
            offset = self.offset_var.get()
            combination = (h, k, l, offset)
            if not self.is_duplicate(combination):
                self.combination_tree.insert('', tk.END, values=combination)
            else:
                messagebox.showwarning("Duplicate Entry", "This hkl and offset already exists.")
        except tk.TclError:
            messagebox.showerror("Invalid input", "Please enter valid integers for h, k, l and a float for offset.")

    def is_duplicate(self, combination):
        combination_converted = tuple((str(i) for i in combination))
        for row_id in self.combination_tree.get_children():
            row_values = self.combination_tree.item(row_id, 'values')
            if tuple(row_values) == combination_converted:
                return True
        return False

    def delete_combination(self):
        selected_item = self.combination_tree.selection()
        if selected_item:
            self.combination_tree.delete(selected_item)
        else:
            messagebox.showwarning("No selection", "Please select a surface to delete.")

    def reset_fields(self):
        self.h_var.set(0)
        self.k_var.set(0)
        self.l_var.set(0)
        self.offset_var.set(0.0)
        self.string_var.set("")
        self.file_var.set("")
        self.browse_button.config(state='normal')

    def calculate(self):
        string_input = self.string_var.get().strip()
        file_input = self.file_var.get().strip()
        if not (string_input or file_input):
            tk.messagebox.showerror("Input Error", "Please provide a refcode or select a file.")
            return

        if not self.combination_tree.get_children():
            tk.messagebox.showerror("Selection Error", "There must be at least one surface in the list.")
            return

        items = self.combination_tree.get_children()
        data = []
        for item in items:
            values = self.combination_tree.item(item, 'values')
            try:
                h = int(values[0])
                k = int(values[1])
                l = int(values[2])
                offset = float(values[3])
                data.append((h, k, l, offset))
            except ValueError as e:
                print(f"Error converting data: {e}")
                continue
        if string_input:
            input_string = string_input  # Use string input if available
        elif file_input:
            input_string = file_input

        output_dir = self.dir_var.get()

        surface_charge_controller = SurfaceChargeController(structure=input_string, output_directory=output_dir,
                                                            hkl_and_offsets=data)
        surface_charge_controller.calculate_surface_charge()
        surface_charge_controller.make_report()
        self.root.destroy()


if __name__ == "__main__":
    if len(sys.argv) > 3 and sys.argv[3].endswith(".m2a"):
        mercury = ApplicationInterface()
        run_from_mercury = True
        input_structure = mercury.input_mol2_file
        app = SurfaceChargeGUI(initial_file_path=input_structure)
        app.root.mainloop()
    else:
        app = SurfaceChargeGUI()
        app.root.mainloop()
