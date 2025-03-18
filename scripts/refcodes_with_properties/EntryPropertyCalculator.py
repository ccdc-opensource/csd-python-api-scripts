#
# This script can be used for any purpose without limitation subject to the
# conditions at http://www.ccdc.cam.ac.uk/Community/Pages/Licences/v2.aspx
#
# This permission notice and the following statement of attribution must be
# included in all copies or substantial portions of this script.
#
# 2025-03-14: created by Jason C. Cole, The Cambridge Crystallographic Data Centre


'''
Utility classes for filtering CSD entries based on a property control file
'''

import ccdc.io

_filter_classes = {}


def register(cls):
    ''' Register a filter class to use in the script.
    :param cls: the class to register.
    '''
    if cls.name() in _filter_classes:
        raise ValueError(f"a class with the name {cls.name()} is already registered. Use a different name")

    _filter_classes[cls.name()] = cls


def filter(name):
    return _filter_classes[name]


def helptext():
    ''' Get help text
    '''
    txt = ""
    for name in _filter_classes.keys():
        cls = _filter_classes[name]
        txt = txt + "    %s -> %s," % (name, cls.helptext())
    return txt[:-1]


class _Filter(object):

    @staticmethod
    def name():
        raise NotImplementedError  # override this

    @staticmethod
    def helptext():
        raise NotImplementedError  # override this

    @staticmethod
    def argument_pair():
        raise NotImplementedError  # override this


class _ComparativeFilter(_Filter):
    def __init__(self, args):
        value = False
        if eval(args.strip()) == 1:
            value = True

        self.expected_value = value

    def value(self):
        raise NotImplementedError  # override this

    def __call__(self, theobject):
        value = self.value(theobject)
        return value == self.expected_value


class _RangeFilter(_Filter):
    def __init__(self, args):
        parts = [p.strip() for p in args.split()]
        self.minimum = float(parts[0])
        self.maximum = float(parts[1])

    def value(self):
        raise NotImplementedError  # override this

    def __call__(self, theobject):
        value = self.value(theobject)
        return value >= self.minimum and value <= self.maximum


class AllowedAtomicNumbersFilter(_Filter):
    def __init__(self, args):
        self.allowed_atomic_numbers = [int(atomic_number) for atomic_number in args.strip().split()]

    @staticmethod
    def name():
        return "allowed atomic numbers"

    @staticmethod
    def helptext():
        return "specify a set of atomic numbers (space separated) that the structure can have (and no others)"

    def __call__(self, entry):
        try:
            molecule = entry.crystal.molecule
            return len([x for x in molecule.atoms if x.atomic_number in self.allowed_atomic_numbers]) == len(
                molecule.atoms)
        except TypeError:
            return False


register(AllowedAtomicNumbersFilter)


class MustContainAtomicNumbersFilter(_Filter):
    def __init__(self, args):
        self.must_have_atomic_numbers = [int(atomic_number) for atomic_number in args.strip().split()]

    @staticmethod
    def name():
        return "must have atomic numbers"

    @staticmethod
    def helptext():
        return "specify a set of atomic numbers (space separated) that the structure must have"

    def __call__(self, entry):
        try:
            molecule = entry.crystal.molecule

            contains = {}
            for x in molecule.atoms:
                if not contains.has_key(x.atomic_number):
                    contains[x.atomic_number] = 0
                contains[x.atomic_number] = contains[x.atomic_number] + 1
            for x in self.must_have_atomic_numbers:
                if not contains.has_key(x):
                    return False

            return True
        except:
            return False


register(MustContainAtomicNumbersFilter)


class OrganicFilter(_ComparativeFilter):
    def __init__(self, args):
        super(self.__class__, self).__init__(args)

    @staticmethod
    def name():
        return "organic"

    @staticmethod
    def helptext():
        return "organic entries or not"

    def value(self, entry):
        return entry.is_organic


register(OrganicFilter)


class PolymericFilter(_ComparativeFilter):
    def __init__(self, args):
        super(self.__class__, self).__init__(args)

    @staticmethod
    def name():
        return "polymeric"

    @staticmethod
    def helptext():
        return "polymeric entries or not"

    def value(self, entry):
        return entry.is_polymeric


register(PolymericFilter)


class AllHaveSitesFilter(_ComparativeFilter):
    def __init__(self, args):
        super(self.__class__, self).__init__(args)

    @staticmethod
    def name():
        return "all atoms have sites"

    @staticmethod
    def helptext():
        return "whether all atoms have to have sites"

    def value(self, entry):
        try:
            return entry.crystal.molecule.all_atoms_have_sites
        except:
            return False


register(AllHaveSitesFilter)


class DisorderedFilter(_ComparativeFilter):
    def __init__(self, args):
        super(self.__class__, self).__init__(args)

    @staticmethod
    def name():
        return "disordered"

    @staticmethod
    def helptext():
        return "disordered entries or not"

    def value(self, entry):
        return entry.has_disorder


register(DisorderedFilter)


class AtomicWeightFilter(_RangeFilter):
    def __init__(self, args):
        super(self.__class__, self).__init__(args)

    @staticmethod
    def name():
        return "atomic weight"

    @staticmethod
    def helptext():
        return "specify a range of atomic weight (for the whole structure - not individual molecules)"

    def value(self, entry):
        try:
            molecule = entry.crystal.molecule
            return molecule.molecular_weight
        except TypeError:
            return 0.0


register(AtomicWeightFilter)


class AtomCountFilter(_RangeFilter):
    def __init__(self, args):
        super(self.__class__, self).__init__(args)

    @staticmethod
    def name():
        return "atom count"

    @staticmethod
    def helptext():
        return "specify a range of atom counts (for the whole structure - not individual molecules)"

    def value(self, entry):
        try:
            molecule = entry.crystal.molecule
            return len(molecule.atoms)
        except TypeError:
            return 0


register(AtomCountFilter)


class RotatableBondFilter(_RangeFilter):
    def __init__(self, args):
        super(self.__class__, self).__init__(args)

    @staticmethod
    def name():
        return "rotatable bond count"

    @staticmethod
    def helptext():
        return "specify the number of rotatable bonds (for the whole structure - not individual molecules)"

    def value(self, entry):
        try:
            molecule = entry.crystal.molecule
            return sum(x.is_rotatable for x in molecule.bonds)
        except TypeError:
            return 0


register(RotatableBondFilter)


class DonorCountFilter(_RangeFilter):
    def __init__(self, args):
        super(self.__class__, self).__init__(args)

    @staticmethod
    def name():
        return "donor count"

    @staticmethod
    def helptext():
        return "specify a donor atom count range (for the whole structure - not individual molecules)"

    def value(self, entry):
        try:
            molecule = entry.crystal.molecule
            return len([x for x in molecule.atoms if x.is_donor])
        except TypeError:
            return 0


register(DonorCountFilter)


class AcceptorCountFilter(_RangeFilter):
    def __init__(self, args):
        super(self.__class__, self).__init__(args)

    @staticmethod
    def name():
        return "acceptor count"

    @staticmethod
    def helptext():
        return "specify an acceptor atom count range (for the whole structure - not individual molecules)"

    def value(self, entry):
        try:
            molecule = entry.crystal.molecule
            return len([x for x in molecule.atoms if x.is_acceptor])
        except TypeError:
            return 0


register(AcceptorCountFilter)


class ComponentCountFilter(_RangeFilter):
    def __init__(self, args):
        super(self.__class__, self).__init__(args)

    @staticmethod
    def name():
        return "component range"

    @staticmethod
    def helptext():
        return "specify a component count range for the whole structure"

    def value(self, entry):
        try:
            return len(entry.crystal.molecule.components)
        except TypeError:
            return 0


register(ComponentCountFilter)


class ZPrimeFilter(_RangeFilter):
    def __init__(self, args):
        super(self.__class__, self).__init__(args)

    @staticmethod
    def name():
        return "zprime range"

    @staticmethod
    def helptext():
        return "specify a z-prime range"

    def value(self, entry):
        return entry.crystal.z_prime


register(ZPrimeFilter)


class RfactorFilter(_RangeFilter):
    def __init__(self, args):
        super(self.__class__, self).__init__(args)

    @staticmethod
    def name():
        return "rfactor range"

    @staticmethod
    def helptext():
        return "specify r-factor range (in %%)"

    def value(self, entry):
        return entry.r_factor


register(RfactorFilter)


class SpacegroupNumberFilter(_RangeFilter):
    def __init__(self, args):
        super(self.__class__, self).__init__(args)

    @staticmethod
    def name():
        return "spacegroup number range"

    @staticmethod
    def helptext():
        return "specify spacegroup number range"

    def value(self, entry):
        return entry.crystal.spacegroup_number_and_setting[0]


register(SpacegroupNumberFilter)


class FilterEvaluation(object):
    def __init__(self):
        self._methods = []

    def add_filter(self, method):
        self._methods.append(method)

    def evaluate(self, entry):
        for method in self._methods:
            try:
                if not method(entry):
                    return False
            except TypeError:
                return False

        return True

    def values(self, entry):
        values = {}
        for method in self._methods:
            if hasattr(method, "value"):
                try:
                    values[method.name()] = method.value(entry)
                except NotImplementedError:
                    pass
        return values


def parse_control_file(lines):
    evaluator = FilterEvaluation()
    for line in lines:
        if len(line) > 0 and line[0] != '#':
            parts = line.split(":")
            if len(parts) > 1:
                cls = _filter_classes[parts[0].strip()]
                evaluator.add_filter(cls(parts[1]))
    return evaluator


import unittest


class TestFiltering(unittest.TestCase):

    def setUp(self):

        self.reader = ccdc.io.EntryReader('CSD')
        self.aabhtz = self.reader.entry("AABHTZ")
        self.aacani_ten = self.reader.entry("AACANI10")
        self.aadamc = self.reader.entry("AADAMC")
        self.aadrib = self.reader.entry("AADRIB")
        self.abadis = self.reader.entry("ABADIS")

    def test_organic_filter(self):

        test_file = """
organic : 1
"""
        lines = test_file.split('\n')
        evaluator = parse_control_file(lines)

        self.assertTrue(evaluator.evaluate(self.aabhtz))

        self.assertFalse(evaluator.evaluate(self.aacani_ten))

    def test_component_filter(self):
        test_file = """
component range : 0 1
"""
        lines = test_file.split('\n')
        evaluator = parse_control_file(lines)

        self.assertTrue(evaluator.evaluate(self.aabhtz))

        self.assertFalse(evaluator.evaluate(self.aacani_ten))

    def test_donor_count_filter(self):
        test_file = """
donor count : 2 2
"""
        lines = test_file.split('\n')
        evaluator = parse_control_file(lines)

        self.assertFalse(evaluator.evaluate(self.aabhtz))

        self.assertTrue(evaluator.evaluate(self.aadamc))

        test_file = """
donor count : 0 3
"""
        lines = test_file.split('\n')
        evaluator = parse_control_file(lines)

        self.assertTrue(evaluator.evaluate(self.aabhtz))
        self.assertTrue(evaluator.evaluate(self.aadamc))

    def test_acceptor_count_filter(self):
        test_file = """
acceptor count : 7 7
"""
        lines = test_file.split('\n')
        evaluator = parse_control_file(lines)

        # regards Cl as an acceptor ...
        self.assertTrue(evaluator.evaluate(self.aabhtz))

        self.assertTrue(evaluator.evaluate(self.aacani_ten))

    def test_zprime(self):
        test_file = """
zprime range : 0.99 1.01
"""
        lines = test_file.split('\n')
        evaluator = parse_control_file(lines)
        self.assertTrue(evaluator.evaluate(self.aabhtz))
        self.assertFalse(evaluator.evaluate(self.aadrib))

    def test_atomic_numbers(self):
        test_file = """
allowed atomic numbers : 1 6 7 8
must have atomic numbers : 1 6 7 8
"""
        lines = test_file.split('\n')
        evaluator = parse_control_file(lines)
        self.assertFalse(evaluator.evaluate(self.aabhtz))
        self.assertFalse(evaluator.evaluate(self.aadrib))

        test_file = """
must have atomic numbers : 1 6 7 8
"""
        lines = test_file.split('\n')
        evaluator = parse_control_file(lines)
        self.assertTrue(evaluator.evaluate(self.aabhtz))
        self.assertFalse(evaluator.evaluate(self.aadrib))

    def test_rotatable_bond_count(self):
        test_file = """
rotatable bond count : 0 3
"""
        lines = test_file.split('\n')
        evaluator = parse_control_file(lines)
        self.assertTrue(evaluator.evaluate(self.abadis))

    def test_multiple(self):
        test_file = """

# An example control file
#
#
# only include organic structures as output
organic : 1
# specify a range of donors
donor count     : 0 10
# specify a range of acceptors
acceptor count :    5 5
# rotatable bond count range
rotatable bond count : 3 7
# number of atoms to allow through
atom count : 0 100
# only include structures containing Hydrogen, Carbon, Nitrogen or Oxygen and nothing else
allowed atomic numbers : 1 6 7 8
# only include structures containing all of these elements (i.e.) Hydrogen, Carbon, Nitrogen or Oxygen
must have atomic numbers : 1 6 7 8
# Ensure Z-prime is one
zprime range : 0.99 1.01
# Ensure only one component in the structure
component range : 2 2
# Dont include disordered structures
disordered : 0
# Specify an R-factor range
rfactor range : 0.1 5
#


"""

        lines = test_file.split('\n')
        evaluator = parse_control_file(lines)

        counter = 0
        hits = []

        test_entries = ['AABHTZ', 'ABAQEB', 'ABELEY', 'ADAQOM', 'ADARAA', 'ADARAZ', 'ADUWIG', 'AFEREK']
        for id in test_entries:
            e = self.reader.entry(id)

            if evaluator.evaluate(e):
                hits.append(e.identifier)

        self.assertEquals(['ABAQEB', 'ABELEY', 'ADAQOM', 'ADUWIG', 'AFEREK'], hits)


if __name__ == "__main__":
    unittest.main()
