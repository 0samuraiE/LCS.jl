#!/bin/bash

# Install Claude Code via official installer (devcontainer features skipped: outdated)
curl -fsSL https://claude.ai/install.sh | bash

# Configure env vars in .zshrc.local to avoid grml-zsh overrides
echo 'export HISTFILE="/root/.zsh_history/history"' >> /root/.zshrc.local
echo 'export PATH="$HOME/.local/bin:$PATH"' >> ~/.zshrc.local

julia --project -e "using Pkg; Pkg.instantiate(); Pkg.precompile()"
julia --project -e "using Schema; Schema.generate()"
