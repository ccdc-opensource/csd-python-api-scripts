import unittest

from ccdc.io import EntryReader

from entry_property_calculator import parse_control_file


class TestFiltering(unittest.TestCase):

    def setUp(self):

        self.reader = EntryReader('CSD')
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
rotatable bond count : 0 4
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

        self.assertEqual(['ABAQEB', 'ABELEY', 'ADAQOM', 'ADUWIG', 'AFEREK'], hits)
