from . import tests 
from . import test_comparison_table
from . import inference_basics

test_index = tests.tests
comparison_table = test_comparison_table.comparison_table
inference_summary = inference_basics.inferential_stats_summary
tailed_tests_differences = inference_basics.tailed_tests_differences
tailed_tests_separate_examples = inference_basics.tailed_tests_separate_examples
tailed_tests_single_example = inference_basics.tailed_tests_single_example

class Catalogue:

    def list_tests(parametric_only=False, non_parametric_only=False):
        if parametric_only:
            return [test for test in test_index.keys() if test_index[test].parametric]
        elif non_parametric_only:
            return [test for test in test_index.keys() if not test_index[test].parametric]
        return test_index.keys()

    def show_usecases():
        return comparison_table

    def get_test(test_name):
        return test_index[test_name]

    def get_test_examples(test_name):
        return test_index[test_name].get('examples')

    def get_test_description(test_name):
        return test_index[test_name].get('description')

    def get_test_formulas(test_name):
        return test_index[test_name].get('formulas')

    def get_test_use_cases(test_name):
        return test_index[test_name].get('use_cases')

    def get_test_summary(test_name):
        return test_index[test_name].get('summary')

    def get_code_snippets(test_name):
        return test_index[test_name].get('code_snippets')

    def get_thorough_examples(test_name):
        return test_index[test_name].get('thorough_examples')

class BasicInference:

    def get_summary():
        return inference_summary

    def get_tailed_tests_difference():
        return tailed_tests_differences

    def get_tailed_tests_separate_examples():
        return tailed_tests_separate_examples

    def get_tailed_tests_single_example():
        return tailed_tests_single_example
