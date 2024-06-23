test: test-unit

test-unit:
	pytest --cov=microdeco --cov-report=xml --cov-report=term

gen-requirements:
	pip freeze > requirements.txt