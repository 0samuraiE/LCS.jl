#!/bin/sh
set -eu

export DEBIAN_FRONTEND=noninteractive

apt-get update
apt-get install -y --no-install-recommends \
  zsh \
  git \
  curl \
  ca-certificates \
  passwd \
  procps \
  util-linux \
  && rm -rf /var/lib/apt/lists/*


ENABLE_COMPLETIONS="${ENABLECOMPLETIONS:-true}"
ENABLE_SYNTAX_HIGHLIGHTING="${ENABLESYNTAXHIGHLIGHTING:-true}"
ENABLE_AUTOSUGGESTIONS="${ENABLEAUTOSUGGESTIONS:-true}"
ENABLE_HISTORY_SUBSTRING="${ENABLEHISTORYSUBSTRING:-true}"

USERNAME="${_REMOTE_USER:-${_CONTAINER_USER:-typst}}"

USER_HOME="$(getent passwd "$USERNAME" | cut -d: -f6 || true)"
[ -n "${USER_HOME:-}" ] || USER_HOME="/home/$USERNAME"
[ -d "$USER_HOME" ] || { echo "no home: $USER_HOME"; exit 1; }

ZSH_DIR="$USER_HOME/.zsh"
mkdir -p "$ZSH_DIR"

curl -fsSL https://grml.org/console/zshrc -o "$USER_HOME/.zshrc"

clone() {
  repo="$1"
  dir="$2"
  [ -d "$dir" ] || git clone --depth=1 "$repo" "$dir"
}

[ "$ENABLE_COMPLETIONS" = "true" ] && clone https://github.com/zsh-users/zsh-completions.git "$ZSH_DIR/zsh-completions"
[ "$ENABLE_SYNTAX_HIGHLIGHTING" = "true" ] && clone https://github.com/zsh-users/zsh-syntax-highlighting.git "$ZSH_DIR/zsh-syntax-highlighting"
[ "$ENABLE_AUTOSUGGESTIONS" = "true" ] && clone https://github.com/zsh-users/zsh-autosuggestions.git "$ZSH_DIR/zsh-autosuggestions"
[ "$ENABLE_HISTORY_SUBSTRING" = "true" ] && clone https://github.com/zsh-users/zsh-history-substring-search.git "$ZSH_DIR/zsh-history-substring-search"

cat > "$USER_HOME/.zshrc.pre" <<'EOF'
fpath+=($HOME/.zsh/zsh-completions)
EOF

cat > "$USER_HOME/.zshrc.local" <<'EOF'
zstyle ':completion:*' menu select=1
[ -f "$HOME/.zsh/zsh-syntax-highlighting/zsh-syntax-highlighting.zsh" ] && . "$HOME/.zsh/zsh-syntax-highlighting/zsh-syntax-highlighting.zsh"
[ -f "$HOME/.zsh/zsh-autosuggestions/zsh-autosuggestions.zsh" ] && . "$HOME/.zsh/zsh-autosuggestions/zsh-autosuggestions.zsh"
[ -f "$HOME/.zsh/zsh-history-substring-search/zsh-history-substring-search.zsh" ] && . "$HOME/.zsh/zsh-history-substring-search/zsh-history-substring-search.zsh"
autoload -U compinit && compinit
EOF

GROUPNAME="$(id -gn "$USERNAME" 2>/dev/null || echo "$USERNAME")"
chown -R "$USERNAME":"$GROUPNAME" \
  "$ZSH_DIR" "$USER_HOME/.zshrc" "$USER_HOME/.zshrc.pre" "$USER_HOME/.zshrc.local" || true

usermod -s /bin/zsh "$USERNAME" 2>/dev/null || true