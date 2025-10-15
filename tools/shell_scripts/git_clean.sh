# Do we still declare it in the current .gitmodules?
git config -f .gitmodules --get-regexp '^submodule\.' || echo "no .gitmodules entries"

# Do we still have a local submodule config for it?
git config --show-origin --get-regexp '^submodule\.' || echo "no local submodule config"

# Is there a leftover submodule repository?
test -d .git/modules/libs/FullContactV2 && echo "HAS cache at .git/modules/libs/FullContactV2" || echo "no submodule cache"

# Which branches/refs still *contain* that submodule path (as a gitlink)?
for r in $(git for-each-ref --format='%(refname:short)' refs/heads refs/remotes); do
  git ls-tree -r "$r" | grep -q 'libs/FullContactV2' && echo "ref contains gitlink: $r"
done
