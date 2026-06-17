#
# This script can be used for any purpose without limitation subject to the
# conditions at http://www.ccdc.cam.ac.uk/Community/Pages/Licences/v2.aspx
#
# This permission notice and the following statement of attribution must be
# included in all copies or substantial portions of this script.
#
# 2016-12-06 - created by Anthony Reilly, The Cambridge Crystallographic Data Centre based on code from Andrew Maloney, CCDC
# 2021-02-24 - updated by Alex Moldovan, The Cambridge Crystallographic Data Centre
#

import os
from ccdc.utilities import ApplicationInterface


def make_castep_input(crystal, kp_spacing=0.06, cut_off=750.000, task='SinglePoint', di_opt=False):
    # cell file
    with open('%s.cell' % crystal.identifier, 'w') as cell_file:
        cell_file.write('%BLOCK LATTICE_ABC\n')
        cell_file.write('  %11.7f  %11.7f  %11.7f\n' % crystal.cell_lengths)
        cell_file.write('  %11.7f  %11.7f  %11.7f\n' % crystal.cell_angles)
        cell_file.write('%ENDBLOCK LATTICE_ABC\n\n')
        cell_file.write('%BLOCK POSITIONS_FRAC\n')
        for op in crystal.symmetry_operators:
            mol = crystal.symmetric_molecule(op, [0, 0, 0], force=False)
            for atom in mol.atoms:
                cell_file.write('  %s  %19.16f  %19.16f  %19.16f\n' % (atom.atomic_symbol,
                                                                       atom.fractional_coordinates.x,
                                                                       atom.fractional_coordinates.y,
                                                                       atom.fractional_coordinates.z))
        cell_file.write('%ENDBLOCK POSITIONS_FRAC\n\n')

        cell_file.write('KPOINT_MP_SPACING ' + str(kp_spacing) + '\n\n')

        cell_file.write('%BLOCK SYMMETRY_OPS\n')
        for op in crystal.symmetry_operators:
            rotation = crystal.symmetry_rotation(op)
            trans = crystal.symmetry_translation(op)
            cell_file.write('      %18.15f      %18.15f      %18.15f\n' % (rotation[0], rotation[1], rotation[2]))
            cell_file.write('      %18.15f      %18.15f      %18.15f\n' % (rotation[3], rotation[4], rotation[5]))
            cell_file.write('      %18.15f      %18.15f      %18.15f\n' % (rotation[6], rotation[7], rotation[8]))
            cell_file.write('      %18.15f      %18.15f      %18.15f\n' % (trans[0], trans[1], trans[2]))
            cell_file.write('###\n')
        cell_file.write('%ENDBLOCK SYMMETRY_OPS\n\n')

    # param file
    with open('%s.param' % crystal.identifier, 'w') as param_file:
        param_file.write('task : ' + task + '\n')
        param_file.write('comment : CASTEP calculation for %s\n' % crystal.identifier)

        param_file.write('xc_functional : PBE\n')
        param_file.write('sedc_scheme : TS\n')

        param_file.write('metals_method : dm\n')
        param_file.write('mixing_scheme : Pulay\n')
        param_file.write('spin_polarized : false\n')
        param_file.write('opt_strategy : speed\n')
        param_file.write('cut_off_energy :  ' + str(cut_off) + '\n')
        param_file.write('grid_scale :        2.0\n')
        param_file.write('fine_grid_scale :   3.0\n')
        param_file.write('elec_energy_tol :   1.000e-008\n')
        param_file.write('fix_occupancy : true\n')
        if task != 'SinglePoint':
            param_file.write('finite_basis_corr : 2\n')

        param_file.write('geom_modulus_est : 50 GPa\n')
        param_file.write('geom_max_iter : 200\n')
        param_file.write('num_backup_iter : 1\n')

        param_file.write('geom_energy_tol : 5E-06\n')
        param_file.write('geom_stress_tol : 0.02\n')
        param_file.write('geom_disp_tol : 1E-03\n')
        param_file.write('geom_force_tol : 5E-03\n')
        if di_opt:
            param_file.write('geom_method : delocalised \n')

        param_file.write('write_cif_structure : true\n')
        param_file.write('write_cell_structure : true\n')
        param_file.write('#continuation : default\n')

    return True


if __name__ == '__main__':
    helper = ApplicationInterface()
    entry = helper.current_entry
    crystal = entry.crystal
    os.chdir(helper.output_directory_path)
    make_castep_input(crystal)
