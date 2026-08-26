#include <dirent.h>
#include <errno.h>
#include <stdio.h>
#include <string.h>

static int has_suffix(const char *value, const char *suffix) {
  size_t value_len = strlen(value);
  size_t suffix_len = strlen(suffix);

  return value_len >= suffix_len &&
         strcmp(value + value_len - suffix_len, suffix) == 0;
}

static int is_mpi_prf_payload(const char *name) {
  if (!has_suffix(name, ".prf")) return 0;
  if (has_suffix(name, "_key.prf")) return 1;
  if (has_suffix(name, "_key_idx.prf")) return 1;
  return strstr(name, "_comp") != NULL && strstr(name, "_chunk_") != NULL;
}

int ff_checkpoint_remove_file(const char *path) {
  if (remove(path) == 0 || errno == ENOENT) return 0;
  return errno;
}

int ff_checkpoint_rename_file(const char *old_path, const char *new_path) {
  if (rename(old_path, new_path) == 0) return 0;
  return errno;
}

int ff_checkpoint_remove_payloads(const char *directory, const char *prefix) {
  DIR *dir = opendir(directory);
  struct dirent *entry;
  size_t prefix_len = strlen(prefix);
  int first_error = 0;

  if (dir == NULL) return errno == ENOENT ? 0 : errno;

  while ((entry = readdir(dir)) != NULL) {
    char path[4096];
    int status;

    if (!is_mpi_prf_payload(entry->d_name)) continue;
    if (prefix_len > 0 &&
        (strncmp(entry->d_name, prefix, prefix_len) != 0 ||
         entry->d_name[prefix_len] != '_')) continue;
    if (snprintf(path, sizeof(path), "%s/%s", directory, entry->d_name) >=
        (int)sizeof(path)) {
      if (first_error == 0) first_error = ENAMETOOLONG;
      continue;
    }
    status = remove(path);
    if (status != 0 && errno != ENOENT && first_error == 0) first_error = errno;
  }

  closedir(dir);
  return first_error;
}

int ff_checkpoint_has_mpi_prf_payload(const char *directory) {
  DIR *dir = opendir(directory);
  struct dirent *entry;
  int found = 0;

  if (dir == NULL) return 0;
  while ((entry = readdir(dir)) != NULL) {
    if (is_mpi_prf_payload(entry->d_name)) {
      found = 1;
      break;
    }
  }
  closedir(dir);
  return found;
}
