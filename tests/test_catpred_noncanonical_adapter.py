from pathlib import Path
import importlib.util
import unittest

spec=importlib.util.spec_from_file_location('adapter',Path(__file__).resolve().parents[1]/'scripts/catpred_noncanonical_entry.py')
adapter=importlib.util.module_from_spec(spec);spec.loader.exec_module(adapter)

class AdapterTests(unittest.TestCase):
    def test_canonical_identity_and_unknown_policy(self):
        source="def tokenize(seq):\n    letter_to_num = dict(zip('ACDEFGHIKLMNPQRSTVWY', range(20)))\n    return [letter_to_num[a] for a in seq]\n"
        before={};after={}
        exec(source,before);exec(adapter.patch_source(source),after)
        standard='ACDEFGHIKLMNPQRSTVWY'
        self.assertEqual(before['tokenize'](standard),after['tokenize'](standard))
        self.assertEqual(after['tokenize']('XUBO'),[20,20,20,20])
        self.assertEqual(len(after['tokenize']('AXC')),3)

    def test_unexpected_upstream_fails_closed(self):
        with self.assertRaises(ValueError):adapter.patch_source('unrelated source')

if __name__=='__main__':unittest.main()
