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
        if args.strip() == '1' or args.strip().lower() == 'true':
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
            return len([x for x in molecule.atoms if x.atomic_number in self.allowed_atomic_numbers]) == len(molecule.atoms)
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

            contains = set()
            for x in molecule.atoms:
                contains.add(x.atomic_number)

            for x in self.must_have_atomic_numbers:
                if x not in contains:
                    return False

            return True
        except TypeError:
            return False


register(MustContainAtomicNumbersFilter)


class OrganicFilter(_ComparativeFilter):
    def __init__(self, args):
        super().__init__(args)

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
        super().__init__(args)

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
        super().__init__(args)

    @staticmethod
    def name():
        return "all atoms have sites"

    @staticmethod
    def helptext():
        return "whether all atoms have to have sites"

    def value(self, entry):
        try:
            return entry.crystal.molecule.all_atoms_have_sites
        except TypeError:
            return False


register(AllHaveSitesFilter)


class DisorderedFilter(_ComparativeFilter):
    def __init__(self, args):
        super().__init__(args)

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
        super().__init__(args)

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
        super().__init__(args)

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
        super().__init__(args)

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
        super().__init__(args)

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
        super().__init__(args)

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
        super().__init__(args)

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
        super().__init__(args)

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
        super().__init__(args)

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
        super().__init__(args)

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
