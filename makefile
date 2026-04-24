help: ## show this help
	@sed -ne "s/^##\(.*\)/\1/p" $(MAKEFILE_LIST)
	@printf "────────────────────────`tput bold``tput setaf 2` Make Commands `tput sgr0`────────────────────────────────\n"
	@sed -ne "/@sed/!s/\(^[^#?=]*:\).*##\(.*\)/`tput setaf 2``tput bold`\1`tput sgr0`\2/p" $(MAKEFILE_LIST)
	@printf "────────────────────────`tput bold``tput setaf 4` Make Variables `tput sgr0`───────────────────────────────\n"
	@sed -ne "/@sed/!s/\(.*\)?=\(.*\)##\(.*\)/`tput setaf 4``tput bold`\1:`tput setaf 5`\2`tput sgr0`\3/p" $(MAKEFILE_LIST)
	@printf "───────────────────────────────────────────────────────────────────────\n"

# remove if you want another default.
.DEFAULT_GOAL := help


start-notebook: ## Start the pluto notebook
	julia -e 'using Pluto;Pluto.run(notebook="notebook.jl", auto_reload_from_file=true)'
