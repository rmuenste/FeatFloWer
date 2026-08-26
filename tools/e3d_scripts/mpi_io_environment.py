"""OpenMPI/ROMIO environment setup shared by the E3D launchers."""

import os
import re
import shutil
import subprocess


ROMIO_ENV_DEFAULTS = {
    "ROMIO_CB_BUFFER_SIZE": "16777216",
    "ROMIO_DS_WRITE": "enable",
}


def detect_openmpi_io_components(which=shutil.which,
                                  check_output=subprocess.check_output):
    """Return ``(status, components)`` for the active OpenMPI installation."""
    ompi_info = which("ompi_info")
    if ompi_info is None:
        return "no_openmpi", ()

    try:
        output = check_output(
            [ompi_info, "--param", "io", "all"],
            stderr=subprocess.STDOUT,
            universal_newlines=True,
        )
    except (OSError, subprocess.CalledProcessError):
        return "query_failed", ()

    components = []
    for component in re.findall(r"MCA io:\s+([A-Za-z0-9_]+)\b", output):
        if component not in components:
            components.append(component)
    if not components:
        return "no_components", ()
    return "available", tuple(components)


def configure_mpi_io_environment(environment=None, emit=print,
                                 detector=detect_openmpi_io_components):
    """Configure safe ROMIO defaults without replacing explicit overrides.

    ``OMPI_MCA_io`` is selected only when it is absent. An explicit component
    must exist in the active OpenMPI installation; this catches stale values
    such as ``romio321`` after an OpenMPI upgrade instead of silently running
    with a different I/O backend.
    """
    environment = os.environ if environment is None else environment
    detection_status, components = detector()

    if detection_status == "no_openmpi":
        emit("OpenMPI not detected (ompi_info not found in PATH); "
             "leaving MPI I/O environment unchanged.")
        return detection_status
    if detection_status == "query_failed":
        emit("OpenMPI detected, but querying its MPI I/O components failed; "
             "leaving MPI I/O environment unchanged.")
        return detection_status
    if detection_status == "no_components":
        emit("OpenMPI detected, but no MPI I/O component is available; "
             "leaving MPI I/O environment unchanged.")
        return detection_status

    configured_component = environment.get("OMPI_MCA_io")
    if configured_component:
        if configured_component not in components:
            raise RuntimeError(
                "OMPI_MCA_io={} is not provided by this OpenMPI installation "
                "(available: {}).".format(
                    configured_component, ", ".join(components)))
        component_status = "kept existing"
    else:
        configured_component = next(
            (name for name in components if name.startswith("romio")), None)
        if configured_component is None:
            emit("OpenMPI detected, but no ROMIO component is available; "
                 "leaving MPI I/O environment unchanged.")
            return "no_romio"
        environment["OMPI_MCA_io"] = configured_component
        component_status = "set detected"

    emit("Configuring MPI environment variables for ROMIO:")
    emit("  Detected OpenMPI I/O components: {}".format(
        ", ".join(components)))
    emit("  OMPI_MCA_io={} [{}]".format(
        configured_component, component_status))

    if not configured_component.startswith("romio"):
        emit("  Selected component is not ROMIO; ROMIO tuning variables "
             "were left unchanged.")
        return "configured_non_romio"

    for variable, default_value in ROMIO_ENV_DEFAULTS.items():
        existing = environment.get(variable)
        if existing is None:
            environment[variable] = default_value
            emit("  {}={} [set default]".format(variable, default_value))
        else:
            emit("  {}={} [kept existing]".format(variable, existing))
    return "configured_romio"
