# Remove from .gitmodules (if present)
git config -f .gitmodules --name-only --get-regexp '^submodule\.' | grep 'FullContactV2' && \
  git config -f .gitmodules --remove-section submodule.libs/FullContactV2 || true

# Stage & commit the cleaned .gitmodules so teammates stop seeing it too
git add .gitmodules 2>/dev/null || true
git commit -m "Cleanup: drop obsolete submodule libs/FullContactV2 from .gitmodules" || true

# Remove any local submodule config section (quotes/path style may differ; list to see exact key)
git config --name-only --get-regexp '^submodule\.' | grep 'FullContactV2' && \
  git config --remove-section submodule.libs/FullContactV2 || true

# Delete leftover cached repo, if it exists
rm -rf .git/modules/libs/FullContactV2