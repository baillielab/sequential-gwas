WORKFLOW = workflow.wdl
PROJECT_ID = project-J0pkqyQJpYQ133JG1p2J1qzv
DXCOMPILER = java -jar ~/dxCompiler/dxCompiler-2.10.2.jar
CROMWELL = java -jar ~/cromwell/cromwell-82.jar

check:
	miniwdl check $(WORKFLOW)

local:
	$(CROMWELL) run --inputs inputs.json $(WORKFLOW)

compile:
	$(DXCOMPILER) compile $(WORKFLOW) \
		-archive \
		-folder /workflows \
		-project $(PROJECT_ID)

clean:
	rm -rf cromwell-workflow-logs cromwell-executions