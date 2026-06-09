#
# This script can be used for any purpose without limitation subject to the
# conditions at https://www.ccdc.cam.ac.uk/Community/Pages/Licences/v2.aspx
#
# This permission notice and the following statement of attribution must be
# included in all copies or substantial portions of this script.
#
# 2025-03-14: created by Jason C. Cole, The Cambridge Crystallographic Data Centre

from abc import ABC, abstractmethod

# Utility classes for filtering CSD entries based on a property control file
_filter_classes = {}


def get_filter(name):
    """Look up a registered filter class by name."""
    return _filter_classes[name]


# Keep backward compatibility
filter = get_filter


def helptext():
    """Get help text"""
    return ", ".join(
        f"    {name} -> {cls.helptext()}" for name, cls in _filter_classes.items()
    )


class _Filter(ABC):
    """Base class for all filters. Subclasses are auto-registered unless their name starts with '_'."""

    def __init_subclass__(cls, **kwargs):
        super().__init_subclass__(**kwargs)
        # Auto-register concrete filters (skip private/abstract base classes)
        if not cls.__name__.startswith('_'):
            try:
                name = cls.name()
                if name in _filter_classes:
                    raise ValueError(
                        f"a class with the name {name} is already registered. Use a different name"
                    )
                _filter_classes[name] = cls
            except (NotImplementedError, TypeError):
                pass  # Abstract intermediate classes

    @staticmethod
    @abstractmethod
    def name():
        ...

    @staticmethod
    @abstractmethod
    def helptext():
        ...

    @abstractmethod
    def __call__(self, entry):
        ...


class _ComparativeFilter(_Filter):
    def __init__(self, args):
        self.expected_value = args.strip() in ('1', 'true', 'True')

    @abstractmethod
    def value(self, theobject):
        ...

    def __call__(self, theobject):
        return self.value(theobject) == self.expected_value


class _RangeFilter(_Filter):
    def __init__(self, args):
        parts = args.split()
        self.minimum = float(parts[0])
        self.maximum = float(parts[1])

    @abstractmethod
    def value(self, theobject):
        ...

    def __call__(self, theobject):
        value = self.value(theobject)
        return self.minimum <= value <= self.maximum


class _ValueFilter(_Filter):
    def __init__(self, args):
        value = args.split()[0]
        self.expected_value = None if value == 'None' else value

    @abstractmethod
    def value(self, theobject):
        ...

    def __call__(self, theobject):
        return self.value(theobject) == self.expected_value


class AllowedAtomicNumbersFilter(_Filter):
    def __init__(self, args):
        self.allowed_atomic_numbers = set(int(n) for n in args.split())

    @staticmethod
    def name():
        return "allowed atomic numbers"

    @staticmethod
    def helptext():
        return "specify a set of atomic numbers (space separated) that the structure can have (and no others)"

    def __call__(self, entry):
        try:
            atoms = entry.crystal.molecule.atoms
            return all(a.atomic_number in self.allowed_atomic_numbers for a in atoms)
        except TypeError:
            return False


class MustContainAtomicNumbersFilter(_Filter):
    def __init__(self, args):
        self.must_have_atomic_numbers = set(int(n) for n in args.split())

    @staticmethod
    def name():
        return "must have atomic numbers"

    @staticmethod
    def helptext():
        return "specify a set of atomic numbers (space separated) that the structure must have"

    def __call__(self, entry):
        try:
            present = {a.atomic_number for a in entry.crystal.molecule.atoms}
            return self.must_have_atomic_numbers <= present
        except TypeError:
            return False


class OrganicFilter(_ComparativeFilter):
    @staticmethod
    def name():
        return "organic"

    @staticmethod
    def helptext():
        return "organic entries or not"

    def value(self, entry):
        return entry.is_organic


class PolymericFilter(_ComparativeFilter):
    @staticmethod
    def name():
        return "polymeric"

    @staticmethod
    def helptext():
        return "polymeric entries or not"

    def value(self, entry):
        return entry.is_polymeric


class AllHaveSitesFilter(_ComparativeFilter):
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


class Has3DStructure(_ComparativeFilter):
    @staticmethod
    def name():
        return "has 3D structure"

    @staticmethod
    def helptext():
        return "whether 3D coordinates have been determined for the structure"

    def value(self, entry):
        return entry.has_3d_structure


class DisorderedFilter(_ComparativeFilter):
    @staticmethod
    def name():
        return "disordered"

    @staticmethod
    def helptext():
        return "disordered entries or not"

    def value(self, entry):
        return entry.has_disorder


class AtomicWeightFilter(_RangeFilter):
    @staticmethod
    def name():
        return "atomic weight"

    @staticmethod
    def helptext():
        return "specify a range of atomic weight (for the whole structure - not individual molecules)"

    def value(self, entry):
        try:
            return entry.crystal.molecule.molecular_weight
        except TypeError:
            return 0.0


class AtomCountFilter(_RangeFilter):
    @staticmethod
    def name():
        return "atom count"

    @staticmethod
    def helptext():
        return "specify a range of atom counts (for the whole structure - not individual molecules)"

    def value(self, entry):
        try:
            return len(entry.crystal.molecule.atoms)
        except TypeError:
            return 0


class RotatableBondFilter(_RangeFilter):
    @staticmethod
    def name():
        return "rotatable bond count"

    @staticmethod
    def helptext():
        return "specify the number of rotatable bonds (for the whole structure - not individual molecules)"

    def value(self, entry):
        try:
            return sum(b.is_rotatable for b in entry.crystal.molecule.bonds)
        except TypeError:
            return 0


class DonorCountFilter(_RangeFilter):
    @staticmethod
    def name():
        return "donor count"

    @staticmethod
    def helptext():
        return "specify a donor atom count range (for the whole structure - not individual molecules)"

    def value(self, entry):
        try:
            return sum(a.is_donor for a in entry.crystal.molecule.atoms)
        except TypeError:
            return 0


class AcceptorCountFilter(_RangeFilter):
    @staticmethod
    def name():
        return "acceptor count"

    @staticmethod
    def helptext():
        return "specify an acceptor atom count range (for the whole structure - not individual molecules)"

    def value(self, entry):
        try:
            return sum(a.is_acceptor for a in entry.crystal.molecule.atoms)
        except TypeError:
            return 0


class ComponentCountFilter(_RangeFilter):
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


class ZPrimeFilter(_RangeFilter):
    @staticmethod
    def name():
        return "zprime range"

    @staticmethod
    def helptext():
        return "specify a z-prime range"

    def value(self, entry):
        return entry.crystal.z_prime


class AsymmUnitFilter(_RangeFilter):
    @staticmethod
    def name():
        return "asymmetric unit components"

    @staticmethod
    def helptext():
        return "specify range of components in the asymmetric unit"

    def value(self, entry):
        return len(entry.crystal.asymmetric_unit_molecule.components)


class RfactorFilter(_RangeFilter):
    @staticmethod
    def name():
        return "rfactor range"

    @staticmethod
    def helptext():
        return "specify r-factor range (in %%)"

    def value(self, entry):
        return entry.r_factor


class SpacegroupNumberFilter(_RangeFilter):
    @staticmethod
    def name():
        return "spacegroup number range"

    @staticmethod
    def helptext():
        return "specify spacegroup number range"

    def value(self, entry):
        return entry.crystal.spacegroup_number_and_setting[0]


class ChiralityFilter(_ValueFilter):
    @staticmethod
    def name():
        return "chirality"

    @staticmethod
    def helptext():
        return "specify the chirality value to be used as filter"

    def value(self, entry):
        try:
            return next(
                (atom.chirality for atom in entry.crystal.molecule.atoms if atom.is_chiral),
                None
            )
        except TypeError:
            return None


class FilterEvaluation:
    def __init__(self):
        self._methods = []

    def add_filter(self, method):
        self._methods.append(method)

    def evaluate(self, entry):
        return all(self._safe_call(m, entry) for m in self._methods)

    @staticmethod
    def _safe_call(method, entry):
        try:
            return method(entry)
        except (TypeError, RuntimeError):
            return False

    def values(self, entry):
        return {
            method.name(): method.value(entry)
            for method in self._methods
            if hasattr(method, "value")
        }


def parse_control_file(lines):
    evaluator = FilterEvaluation()
    for line in lines:
        if line and not line.startswith('#'):
            filter_name, sep, args = line.partition(":")
            if sep:
                cls = _filter_classes[filter_name.strip()]
                evaluator.add_filter(cls(args))
    return evaluator
