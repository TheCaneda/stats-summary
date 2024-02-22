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

    @staticmethod
    def list_tests(parametric_only=False, non_parametric_only=False):
        # Assuming test_index is defined elsewhere and accessible
        if parametric_only:
            return [test for test in test_index.keys() if test_index[test].parametric]
        elif non_parametric_only:
            return [test for test in test_index.keys() if not test_index[test].parametric]
        return list(test_index.keys())

    @staticmethod
    def show_usecases():
        # Assuming comparison_table is defined elsewhere and accessible
        print(comparison_table)

    @staticmethod
    def get_test(test_name):
        # Ensure string conversion or handling for non-string types
        test_details = str(test_index[test_name])
        print(f'[{test_name.upper()}]\n{test_details}')

    @staticmethod
    def get_test_examples(test_name):
        examples = "\n".join(test_index[test_name].get('examples', []))  # Convert list to string
        print(f'[{test_name.upper()} EXAMPLES]\n{examples}')

    @staticmethod
    def get_test_description(test_name):
        description = str(test_index[test_name].get('description', ''))
        print(f'[{test_name.upper()} DESCRIPTION]\n{description}')

    @staticmethod
    def get_test_formulas(test_name):
        formulas = "\n".join(test_index[test_name].get('formulas', []))  # Convert list to string
        print(f'[{test_name.upper()} FORMULAS]\n{formulas}')

    @staticmethod
    def get_test_use_cases(test_name):
        use_cases = "\n".join(test_index[test_name].get('use_cases', []))  # Convert list to string
        print(f'[{test_name.upper()} USE CASES]\n{use_cases}')

    @staticmethod
    def get_test_summary(test_name):
        summary = str(test_index[test_name].get('summary', ''))
        print(f'[{test_name.upper()} SUMMARY]\n{summary}')

    @staticmethod
    def get_code_snippets(test_name):
        code_snippets = "\n".join(test_index[test_name].get('code_snippets', []))  # Convert list to string
        print(f'[{test_name.upper()} CODE SNIPPETS]\n{code_snippets}')

    @staticmethod
    def get_thorough_examples(test_name):
        thorough_examples = "\n".join(test_index[test_name].get('thorough_examples', []))  # Convert list to string
        print(f"[{test_name.upper()} THOROUGH EXAMPLES]\n{thorough_examples}")

class BasicInference:

    def get_summary():
        print(inference_summary)

    def get_tailed_tests_difference():
        print(tailed_tests_differences)

    def get_tailed_tests_separate_examples():
        print(tailed_tests_separate_examples)

    def get_tailed_tests_single_example():
        print(tailed_tests_single_example)
