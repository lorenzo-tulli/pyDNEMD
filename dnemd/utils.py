import logging
import subprocess
import sys
from pathlib import Path


def get_logger(name: str) -> logging.Logger:
    logger = logging.getLogger(name)
    if not logger.handlers:
        h = logging.StreamHandler(sys.stdout)
        h.setFormatter(logging.Formatter("[%(asctime)s %(levelname)s] %(message)s", "%H:%M:%S"))
        logger.addHandler(h)
    logger.setLevel(logging.INFO)
    return logger


def run(cmd: list[str], cwd: str | Path = None, check: bool = True) -> subprocess.CompletedProcess:
    logger = get_logger("run")
    logger.info("$ " + " ".join(str(c) for c in cmd))
    result = subprocess.run(
        [str(c) for c in cmd],
        cwd=str(cwd) if cwd else None,
        check=check,
        text=True,
        capture_output=True,
    )
    if result.stdout:
        logger.info(result.stdout[-2000:])   # tail to avoid flooding logs
    if result.stderr:
        logger.warning(result.stderr[-2000:])
    return result


def run_piped(cmd: list[str], stdin_text: str, cwd: str | Path = None) -> subprocess.CompletedProcess:
    """Run a command with stdin text (for interactive GROMACS prompts)."""
    logger = get_logger("run_piped")
    logger.info("$ " + " ".join(str(c) for c in cmd))
    result = subprocess.run(
        [str(c) for c in cmd],
        input=stdin_text,
        cwd=str(cwd) if cwd else None,
        check=True,
        text=True,
        capture_output=True,
    )
    if result.stdout:
        logger.info(result.stdout[-2000:])
    if result.stderr:
        logger.warning(result.stderr[-2000:])
    return result


def ensure_dir(path: str | Path) -> Path:
    p = Path(path)
    p.mkdir(parents=True, exist_ok=True)
    return p


def copy_file(src: str | Path, dst: str | Path):
    import shutil
    shutil.copy2(str(src), str(dst))
