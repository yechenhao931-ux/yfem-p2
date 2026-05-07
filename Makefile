.PHONY: clean

clean:
	@echo "delete all build/"
	find . -type d -name "build" -exec rm -rf {} +
	@echo "completed"
