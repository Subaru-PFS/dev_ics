import argparse
import logging
import os
import pathlib
import sys
import warnings
from datetime import datetime

import h5py
import numpy as np
import pandas as pd
from tqdm.auto import tqdm

# Load site-specific modules only after the environment's standard packages.
sys.path.append(f'{os.environ["HOME"]}/pythonmodules')
sys.path.append(f'{os.environ["HOME"]}/pythonmodules/selenium')

from ics.cobraCharmer.cobraCoach import visDianosticPlot_TK as vdp


DEFAULT_XML = (
    "/software/mhs/products/Linux64/"
    "pfs_instdata/1.7.36/data/pfi/modules/ALL/ALL.xml"
)


def runDir2Date(runDir):
    """Return the date encoded by a yyyymmdd_numeric run directory."""

    try:
        dateText, runNumber = str(runDir).split("_", maxsplit=1)
        if not runNumber or not runNumber.isdigit():
            raise ValueError
        return datetime.strptime(dateText, "%Y%m%d")
    except ValueError as exc:
        raise ValueError(
            f"Invalid runDir {runDir!r}; expected yyyymmdd followed "
            "by an underscore and a nonempty numeric run number."
        ) from exc


def normalizeEL(value):
    """Return a finite numeric elevation, accepting values such as EL90."""

    if isinstance(value, bytes):
        value = value.decode("utf-8")
    if isinstance(value, str) and value.upper().startswith("EL"):
        value = value[2:]

    try:
        elevation = float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"Invalid EL value: {value!r}") from exc

    if not np.isfinite(elevation):
        raise ValueError(f"Invalid EL value: {value!r}")
    return int(elevation) if elevation.is_integer() else elevation


def askYesNo(prompt):
    while True:
        answer = input(prompt).strip().lower()
        if answer in {"y", "yes"}:
            return True
        if answer in {"n", "no"}:
            return False
        print("Please type yes or no.")


def validateArm(arm):
    arm = str(arm).lower()
    if arm not in {"theta", "phi"}:
        raise ValueError("arm must be either 'theta' or 'phi'.")
    return arm


def visit2runDir(
    visit,
    arm,
    EL,
    *,
    visitRunDirList="visit_runDir.csv",
    update=True,
):
    """Resolve a visit and maintain its runDir, arm, and EL metadata."""

    arm = validateArm(arm)
    EL = normalizeEL(EL)
    csvPath = pathlib.Path(visitRunDirList)

    try:
        table = pd.read_csv(csvPath)
    except FileNotFoundError as exc:
        raise FileNotFoundError(f"Could not find {csvPath}.") from exc

    required = {"visit", "runDir", "arm", "EL"}
    missing = required - set(table.columns)
    if missing:
        raise ValueError(f"{csvPath} is missing columns: {sorted(missing)}")

    visitMask = table["visit"].eq(visit)
    matches = table.loc[visitMask]
    if len(matches) == 1:
        runDir = str(matches.iloc[0]["runDir"])
    else:
        if len(matches) > 1:
            warnings.warn(
                f"Multiple runDir entries exist for visit {visit}; "
                "resolving it again with vdp.findRunDir().",
                RuntimeWarning,
            )
        runDir = vdp.findRunDir(visit)
        if runDir is None:
            raise FileNotFoundError(
                f"No runDir could be found for visit {visit}."
            )
        runDir = str(runDir)

    runDate = runDir2Date(runDir)
    selectedMask = visitMask & table["runDir"].astype(str).eq(runDir)
    changed = False

    if not selectedMask.any():
        newRow = {column: np.nan for column in table.columns}
        newRow.update(visit=visit, runDir=runDir, arm=arm, EL=EL)
        if "year" in table.columns:
            newRow["year"] = runDate.year
        if "month" in table.columns:
            newRow["month"] = runDate.month
        table = pd.concat(
            [table, pd.DataFrame([newRow])], ignore_index=True
        )
        selectedMask = (
            table["visit"].eq(visit)
            & table["runDir"].astype(str).eq(runDir)
        )
        changed = True

    selected = table.loc[selectedMask].iloc[0]
    oldArm = str(selected["arm"]).lower()
    oldEL = normalizeEL(selected["EL"])

    if oldArm != arm and (
        update
        or askYesNo(
            f"CSV arm is {oldArm!r}, requested arm is {arm!r}. "
            "Update it? [y/n] "
        )
    ):
        table.loc[selectedMask, "arm"] = arm
        changed = True

    if oldEL != EL and (
        update
        or askYesNo(
            f"CSV EL is {oldEL!r}, requested EL is {EL!r}. "
            "Update it? [y/n] "
        )
    ):
        table.loc[selectedMask, "EL"] = EL
        changed = True

    duplicateMask = visitMask & ~table["runDir"].astype(str).eq(runDir)
    if duplicateMask.any() and (
        update
        or askYesNo(
            f"Remove other runDir rows for visit {visit}? [y/n] "
        )
    ):
        table = table.loc[~duplicateMask].copy()
        changed = True

    if changed:
        table.to_csv(csvPath, index=False)
    return runDir


def visitList2runDirList(
    visitsFile,
    *,
    visitRunDirList="visit_runDir.csv",
    update=True,
    errorHandler=None,
):
    """Read `runNumber visit arm EL rot` rows and resolve runDirs."""

    table = np.loadtxt(visitsFile, dtype=str, ndmin=2)
    if table.shape[1] != 5:
        raise ValueError(
            f"{visitsFile} must contain: "
            "runNumber, visit, arm, EL, and rot."
        )

    count = len(table)
    runNumbers = np.full(count, None, dtype=object)
    runDirs = np.full(count, None, dtype=object)
    visits = np.full(count, None, dtype=object)
    arms = np.full(count, None, dtype=object)
    elevations = np.full(count, np.nan, dtype=float)
    rotations = np.full(count, np.nan, dtype=float)

    for index, row in enumerate(table):
        try:
            runNumber = int(row[0])
            visit = int(row[1])
            arm = validateArm(row[2])
            EL = normalizeEL(row[3])
            rot = normalizeEL(row[4])
            runDirs[index] = visit2runDir(
                visit,
                arm,
                EL,
                visitRunDirList=visitRunDirList,
                update=update,
            )
            runNumbers[index] = runNumber
            visits[index] = visit
            arms[index] = arm
            elevations[index] = EL
            rotations[index] = rot
        except Exception as exc:
            item = f"row {index + 1} ({' '.join(row)})"
            if errorHandler is None:
                warnings.warn(
                    f"Could not process {item}: {exc}", RuntimeWarning
                )
            else:
                errorHandler(item, exc)

    return runNumbers, runDirs, visits, arms, elevations, rotations


def positionsToAngles(posFW, posRV, center):
    """Convert repeated forward/reverse positions to unwrapped angles."""

    posFW = np.asarray(posFW)
    posRV = np.asarray(posRV)
    if posFW.ndim != 2 or posRV.ndim != 2:
        raise ValueError("Forward and reverse positions must be 2-D.")
    if posFW.shape != posRV.shape:
        raise ValueError("Forward and reverse position shapes do not match.")

    angFW = np.angle(posFW - center).astype(float)
    angRV = np.angle(posRV - center).astype(float)

    for repeatIndex in range(posFW.shape[0]):
        home = angFW[repeatIndex, 0]
        angFW[repeatIndex] = (angFW[repeatIndex] - home) % (2 * np.pi)
        angRV[repeatIndex] = (angRV[repeatIndex] - home) % (2 * np.pi)

        fw = angFW[repeatIndex]
        fwMid = int(np.argmin(np.abs(fw % (2 * np.pi) - np.pi)))
        leading = np.arange(fwMid)
        fw[leading[fw[leading] > 5.0]] -= 2 * np.pi
        if fwMid < fw.size - 1:
            trailing = np.arange(fwMid, fw.size)
            fw[trailing[fw[trailing] < 1.0]] += 2 * np.pi

        rv = angRV[repeatIndex]
        rvMid = int(np.argmin(np.abs(rv % (2 * np.pi) - np.pi)))
        leading = np.arange(rvMid)
        rv[leading[rv[leading] < 1.0]] += 2 * np.pi
        if rvMid < rv.size - 1:
            trailing = np.arange(rvMid, rv.size)
            rv[trailing[rv[trailing] > 5.0]] -= 2 * np.pi
        if rv[0] < 1.0:
            rv += 2 * np.pi

    return angFW, angRV


def calculateMotorMap(
    angFW,
    angRV,
    *,
    binSize=np.deg2rad(3.6),
    regions=112,
    delta=0.1,
    steps=500,
):
    """Calculate forward and reverse radians-per-step motor maps."""

    if steps <= 0:
        raise ValueError("steps must be positive.")
    if binSize <= 0 or regions <= 0:
        raise ValueError("binSize and regions must be positive.")

    angFW = np.asarray(angFW, dtype=float)
    angRV = np.asarray(angRV, dtype=float)
    if angFW.shape != angRV.shape or angFW.ndim != 2:
        raise ValueError("angFW and angRV must be matching 2-D arrays.")

    mmFW = np.zeros(regions, dtype=float)
    mmRV = np.zeros(regions, dtype=float)
    repeatCount, positionCount = angFW.shape

    for binIndex in range(regions):
        binMin = binSize * binIndex
        binMax = binMin + binSize

        fractionSum = 0.0
        valueSum = 0.0
        for repeatIndex in range(repeatCount):
            for positionIndex in range(positionCount - 1):
                start = angFW[repeatIndex, positionIndex]
                end = angFW[repeatIndex, positionIndex + 1]
                if not np.isfinite(start) or not np.isfinite(end):
                    continue
                if end <= start:
                    continue
                if start < binMax and end > binMin:
                    moveInBin = min(end, binMax) - max(start, binMin)
                    entireMove = end - start
                    fraction = moveInBin**2 / entireMove
                    fractionSum += fraction
                    valueSum += fraction * entireMove / steps
        if fractionSum > 0:
            mmFW[binIndex] = valueSum / fractionSum

        fractionSum = 0.0
        valueSum = 0.0
        for repeatIndex in range(repeatCount):
            for positionIndex in range(positionCount - 1):
                start = angRV[repeatIndex, positionIndex]
                end = angRV[repeatIndex, positionIndex + 1]
                if not np.isfinite(start) or not np.isfinite(end):
                    continue
                if end >= start or end < delta:
                    continue
                if start > binMin and end < binMax:
                    moveInBin = min(start, binMax) - max(end, binMin)
                    entireMove = start - end
                    fraction = moveInBin**2 / entireMove
                    fractionSum += fraction
                    valueSum += fraction * entireMove / steps
        if fractionSum > 0:
            mmRV[binIndex] = valueSum / fractionSum

    bad = False
    nonzero = np.flatnonzero(mmFW)
    if nonzero.size:
        mmFW[: nonzero[0]] = mmFW[nonzero[0]]
        mmFW[nonzero[-1] + 1 :] = mmFW[nonzero[-1]]
    else:
        bad = True

    nonzero = np.flatnonzero(mmRV)
    if nonzero.size:
        mmRV[: nonzero[0]] = mmRV[nonzero[0]]
        mmRV[nonzero[-1] + 1 :] = mmRV[nonzero[-1]]
    else:
        bad = True

    return mmFW, mmRV, bad


def loadSteps(runDir, *, dataRoot="/data/MCS", steps=None):
    """Return a caller-supplied step count or load data/steps.npy."""

    if steps is None:
        stepsPath = pathlib.Path(dataRoot) / str(runDir) / "data/steps.npy"
        try:
            loaded = np.asarray(np.load(stepsPath)).squeeze()
        except FileNotFoundError as exc:
            raise FileNotFoundError(f"Could not find {stepsPath}.") from exc
        if loaded.size != 1:
            raise ValueError(f"{stepsPath} must contain one step count.")
        steps = int(loaded)
    else:
        steps = int(steps)

    if steps <= 0:
        raise ValueError(f"Invalid step count for {runDir}: {steps}.")
    return steps


def measureRunMotorMap(
    runDir,
    arm,
    EL,
    *,
    runNumber=None,
    rot=None,
    xml=DEFAULT_XML,
    thetaGeometryHdf5="geometry/thetaGeometry.hdf5",
    phiGeometryHdf5="geometry/phiGeometry.hdf5",
    motorMapHdf5="motorMaps/motorMaps.hdf5",
    dataRoot="/data/MCS",
    binSizeDegrees=3.6,
    regions=112,
    delta=0.1,
    steps=None,
    update=True,
):
    """Measure one arm and save it under motorMaps.hdf5/runDir/arm."""

    runDir = str(runDir)
    runDir2Date(runDir)
    arm = validateArm(arm)
    EL = normalizeEL(EL)
    if (runNumber is None) != (rot is None):
        raise ValueError("runNumber and rot must be supplied together.")
    if runNumber is not None:
        runNumber = int(runNumber)
        rot = normalizeEL(rot)
    geometryFile = pathlib.Path(
        thetaGeometryHdf5 if arm == "theta" else phiGeometryHdf5
    )

    try:
        vis = vdp.VisDianosticPlot(
            runDir=runDir,
            xml=pathlib.Path(xml),
            arm=arm,
            datatype="MM",
        )
        forward = np.asarray(vis.fw)
        reverse = np.asarray(vis.rv)
    except Exception as exc:
        raise AttributeError(f"Cannot import {runDir} Cobra data.") from exc

    if forward.shape != reverse.shape or forward.ndim != 3:
        raise ValueError(
            f"Unexpected position shape for {runDir}: "
            f"fw={forward.shape}, rv={reverse.shape}."
        )

    try:
        with h5py.File(geometryFile, "r") as geometry:
            if runDir not in geometry:
                raise KeyError(
                    f"No geometry group exists for runDir {runDir}."
                )
            group = geometry[runDir]
            missing = {"cobraIdx", "newCnts_circ"} - set(group.keys())
            if missing:
                raise KeyError(
                    f"{geometryFile}/{runDir} is missing {sorted(missing)}."
                )
            cobraIdx = np.asarray(group["cobraIdx"][:])
            centers = np.asarray(group["newCnts_circ"][:], dtype=float)
    except FileNotFoundError as exc:
        raise FileNotFoundError(
            f"Geometry file {geometryFile} does not exist."
        ) from exc

    cobraCount = forward.shape[0]
    if cobraIdx.ndim != 1 or len(cobraIdx) != cobraCount:
        raise ValueError(
            f"cobraIdx in {geometryFile}/{runDir} does not match "
            f"the {cobraCount} measured cobras."
        )
    if centers.shape != (cobraCount, 2):
        raise ValueError(
            f"newCnts_circ in {geometryFile}/{runDir} must have "
            f"shape ({cobraCount}, 2), not {centers.shape}."
        )

    stepCount = loadSteps(runDir, dataRoot=dataRoot, steps=steps)
    binSize = np.deg2rad(float(binSizeDegrees))
    binMin = np.arange(regions, dtype=float) * binSize
    binMax = binMin + binSize
    mmFW = np.full((cobraCount, regions), np.nan)
    mmRV = np.full((cobraCount, regions), np.nan)
    bad = np.ones(cobraCount, dtype=bool)

    previousLoggingLevel = logging.root.manager.disable
    logging.disable(logging.CRITICAL)
    try:
        for cobraRow in tqdm(
            range(cobraCount),
            desc=f"Measuring {arm} motor maps for {runDir}",
        ):
            centerXY = centers[cobraRow]
            if not np.all(np.isfinite(centerXY)):
                continue
            center = centerXY[0] + 1j * centerXY[1]
            try:
                angFW, angRV = positionsToAngles(
                    forward[cobraRow], reverse[cobraRow], center
                )
                cobraFW, cobraRV, cobraBad = calculateMotorMap(
                    angFW,
                    angRV,
                    binSize=binSize,
                    regions=regions,
                    delta=delta,
                    steps=stepCount,
                )
            except Exception as exc:
                warnings.warn(
                    f"Cobra {cobraIdx[cobraRow]} in {runDir} failed: "
                    f"{exc}",
                    RuntimeWarning,
                )
                continue

            mmFW[cobraRow] = cobraFW
            mmRV[cobraRow] = cobraRV
            bad[cobraRow] = cobraBad
    finally:
        logging.disable(previousLoggingLevel)

    outputPath = pathlib.Path(motorMapHdf5)
    outputPath.parent.mkdir(parents=True, exist_ok=True)
    with h5py.File(outputPath, "a") as output:
        runGroup = output.require_group(runDir)
        if arm in runGroup:
            replace = update or askYesNo(
                f"Update {runDir}/{arm} in {outputPath}? [y/n] "
            )
            if not replace:
                return runDir

        if "EL" in runGroup.attrs:
            oldEL = normalizeEL(runGroup.attrs["EL"])
            if oldEL != EL:
                replaceEL = update or askYesNo(
                    f"Stored EL is {oldEL!r}, requested EL is {EL!r}. "
                    "Update it? [y/n] "
                )
                if not replaceEL:
                    return runDir

        runGroup.attrs["runDir"] = runDir
        runGroup.attrs["EL"] = EL
        if runNumber is not None:
            runGroup.attrs["runNumber"] = runNumber
            runGroup.attrs["rot"] = rot
        if arm in runGroup:
            del runGroup[arm]

        armGroup = runGroup.create_group(arm)
        armGroup.attrs["arm"] = arm
        armGroup.attrs["steps"] = stepCount
        armGroup.attrs["binSizeDegrees"] = float(binSizeDegrees)
        armGroup.attrs["regions"] = int(regions)
        armGroup.attrs["delta"] = float(delta)
        armGroup.create_dataset("cobraIdx", data=cobraIdx)
        armGroup.create_dataset("binMin", data=binMin)
        armGroup.create_dataset("binMax", data=binMax)
        armGroup.create_dataset("mmFW", data=mmFW)
        armGroup.create_dataset("mmRV", data=mmRV)
        armGroup.create_dataset("bad", data=bad)

    return runDir


def collectMotorMapsByRunELRotArm(motorMapHdf5):
    """Collect maps into (runNumber, EL, rot, arm) bins."""

    grouped = {}
    with h5py.File(motorMapHdf5, "r") as motorMaps:
        for runDir, runGroup in motorMaps.items():
            if not isinstance(runGroup, h5py.Group):
                continue
            missingAttrs = {
                "runNumber", "EL", "rot"
            } - set(runGroup.attrs.keys())
            if missingAttrs:
                warnings.warn(
                    f"Skipping {motorMapHdf5}/{runDir}: missing attributes "
                    f"{sorted(missingAttrs)}.",
                    RuntimeWarning,
                )
                continue
            runNumber = int(runGroup.attrs["runNumber"])
            EL = normalizeEL(runGroup.attrs["EL"])
            rot = normalizeEL(runGroup.attrs["rot"])

            for arm in ("theta", "phi"):
                if arm not in runGroup:
                    continue
                armGroup = runGroup[arm]
                required = {
                    "cobraIdx", "binMin", "binMax", "mmFW", "mmRV", "bad"
                }
                missing = required - set(armGroup.keys())
                if missing:
                    raise KeyError(
                        f"{motorMapHdf5}/{runDir}/{arm} is missing "
                        f"{sorted(missing)}."
                    )

                cobraIdx = np.asarray(armGroup["cobraIdx"][:])
                mmFW = np.asarray(armGroup["mmFW"][:], dtype=float)
                mmRV = np.asarray(armGroup["mmRV"][:], dtype=float)
                bad = np.asarray(armGroup["bad"][:], dtype=bool)
                if cobraIdx.ndim != 1 or len(np.unique(cobraIdx)) != len(cobraIdx):
                    raise ValueError(
                        f"{motorMapHdf5}/{runDir}/{arm}/cobraIdx is invalid."
                    )
                if mmFW.shape != mmRV.shape or mmFW.shape[0] != len(cobraIdx):
                    raise ValueError(
                        f"Motor-map shapes do not match cobraIdx for "
                        f"{runDir}/{arm}."
                    )
                if bad.shape != (len(cobraIdx),):
                    raise ValueError(
                        f"bad does not match cobraIdx for {runDir}/{arm}."
                    )

                grouped.setdefault((runNumber, EL, rot, arm), []).append(
                    {
                        "runDir": str(runDir),
                        "date": runDir2Date(runDir),
                        "cobraIdx": cobraIdx,
                        "binMin": np.asarray(armGroup["binMin"][:], dtype=float),
                        "binMax": np.asarray(armGroup["binMax"][:], dtype=float),
                        "mmFW": mmFW,
                        "mmRV": mmRV,
                        "bad": bad,
                    }
                )
    return grouped


def medianMotorMaps(runs):
    """Median-combine good run maps, aligned by cobraIdx."""

    referenceMin = runs[0]["binMin"]
    referenceMax = runs[0]["binMax"]
    for run in runs[1:]:
        if not (
            np.array_equal(run["binMin"], referenceMin)
            and np.array_equal(run["binMax"], referenceMax)
        ):
            raise ValueError(
                "Cannot merge motor maps with different bin boundaries."
            )

    cobraIdx = np.unique(
        np.concatenate([run["cobraIdx"] for run in runs])
    )
    regions = len(referenceMin)
    fwStack = np.full((len(runs), len(cobraIdx), regions), np.nan)
    rvStack = np.full_like(fwStack, np.nan)
    validStack = np.zeros((len(runs), len(cobraIdx)), dtype=bool)

    for runIndex, run in enumerate(runs):
        positions = np.searchsorted(cobraIdx, run["cobraIdx"])
        valid = (
            ~run["bad"]
            & np.any(np.isfinite(run["mmFW"]), axis=1)
            & np.any(np.isfinite(run["mmRV"]), axis=1)
        )
        validPositions = positions[valid]
        fwStack[runIndex, validPositions] = run["mmFW"][valid]
        rvStack[runIndex, validPositions] = run["mmRV"][valid]
        validStack[runIndex, validPositions] = True

    with warnings.catch_warnings():
        warnings.filterwarnings(
            "ignore",
            message="All-NaN slice encountered",
            category=RuntimeWarning,
        )
        mmFW = np.nanmedian(fwStack, axis=0)
        mmRV = np.nanmedian(rvStack, axis=0)

    numberOfValidRuns = np.sum(validStack, axis=0)
    bad = numberOfValidRuns == 0
    return {
        "cobraIdx": cobraIdx,
        "binMin": referenceMin,
        "binMax": referenceMax,
        "mmFW": mmFW,
        "mmRV": mmRV,
        "bad": bad,
        "numberOfValidRuns": numberOfValidRuns,
    }


def mergeMotorMapsByRunELRotArm(
    motorMapHdf5,
    combinedMotorMapHdf5,
    *,
    overwrite=True,
):
    """Median-combine maps by run number, EL, rotation, and arm."""

    inputPath = pathlib.Path(motorMapHdf5)
    if not inputPath.exists():
        raise FileNotFoundError(f"Motor-map file {inputPath} does not exist.")
    grouped = collectMotorMapsByRunELRotArm(inputPath)

    outputPath = pathlib.Path(combinedMotorMapHdf5)
    outputPath.parent.mkdir(parents=True, exist_ok=True)
    mode = "w" if overwrite else "x"
    with h5py.File(outputPath, mode) as output:
        runDates = {}
        for (runNumber, _EL, _rot, _arm), runs in grouped.items():
            runDates.setdefault(runNumber, []).extend(
                run["date"] for run in runs
            )

        for runNumber, dates in runDates.items():
            startDate = min(dates)
            endDate = max(dates)
            middleDate = startDate + (endDate - startDate) / 2
            runGroup = output.require_group(f"run{runNumber}")
            runGroup.attrs["runNumber"] = runNumber
            runGroup.attrs["date"] = middleDate.strftime("%Y%m")
            runGroup.attrs["startDate"] = startDate.strftime("%Y%m%d")
            runGroup.attrs["endDate"] = endDate.strftime("%Y%m%d")

        for runNumber, EL, rot, arm in sorted(
            grouped,
            key=lambda key: (
                int(key[0]), float(key[1]), float(key[2]), key[3]
            ),
        ):
            runs = grouped[(runNumber, EL, rot, arm)]
            medians = medianMotorMaps(runs)
            runGroup = output[f"run{runNumber}"]
            elGroup = runGroup.require_group(f"EL_{EL}")
            elGroup.attrs["EL"] = EL
            rotGroup = elGroup.require_group(f"rot_{rot}")
            rotGroup.attrs["rot"] = rot
            armGroup = rotGroup.create_group(arm)
            armGroup.attrs["arm"] = arm
            armGroup.attrs["numberOfRuns"] = len(runs)
            armGroup.attrs["runDirs"] = np.asarray(
                [run["runDir"] for run in runs],
                dtype=h5py.string_dtype("utf-8"),
            )
            for name, values in medians.items():
                armGroup.create_dataset(name, data=values)

    return outputPath


class MotorMapMeasurement:
    def __init__(
        self,
        visit=None,
        arm=None,
        EL=None,
        runNumber=None,
        rot=None,
        runDir=None,
        *,
        xml=DEFAULT_XML,
        visitRunDirList="visit_runDir.csv",
        thetaGeometryHdf5="geometry/thetaGeometry.hdf5",
        phiGeometryHdf5="geometry/phiGeometry.hdf5",
        motorMapHdf5="motorMaps/motorMaps.hdf5",
        combinedMotorMapHdf5="motorMaps/combinedMotorMaps.hdf5",
        dataRoot="/data/MCS",
        binSizeDegrees=3.6,
        regions=112,
        delta=0.1,
        steps=None,
        update=True,
        overwriteCombined=True,
        batchErrorPolicy="warn",
    ):
        if visit is not None and runDir is not None:
            raise ValueError("Specify either visit or runDir, not both.")
        if batchErrorPolicy not in {"warn", "prompt"}:
            raise ValueError(
                "batchErrorPolicy must be either 'warn' or 'prompt'."
            )

        self.visit = visit
        self.arm = str(arm).lower() if arm is not None else None
        self.EL = EL
        self.runNumber = runNumber
        self.rot = rot
        self.runDir = runDir
        self.xml = pathlib.Path(xml)
        self.visitRunDirList = pathlib.Path(visitRunDirList)
        self.thetaGeometryHdf5 = pathlib.Path(thetaGeometryHdf5)
        self.phiGeometryHdf5 = pathlib.Path(phiGeometryHdf5)
        self.motorMapHdf5 = pathlib.Path(motorMapHdf5)
        self.combinedMotorMapHdf5 = pathlib.Path(combinedMotorMapHdf5)
        self.dataRoot = pathlib.Path(dataRoot)
        self.binSizeDegrees = float(binSizeDegrees)
        self.regions = int(regions)
        self.delta = float(delta)
        self.steps = steps
        self.update = update
        self.overwriteCombined = overwriteCombined
        self.batchErrorPolicy = batchErrorPolicy
        self.errors = []
        self.runDirs = None
        self.runNumbers = None
        self.visits = None
        self.arms = None
        self.ELs = None
        self.rots = None

    def _validateSingleInput(self):
        self.arm = validateArm(self.arm)
        if self.EL is None:
            raise ValueError("EL is required for a single visit/runDir.")
        if self.runNumber is None or self.rot is None:
            raise ValueError(
                "runNumber and rot are required for a single visit/runDir."
            )

    def _handleBatchError(self, item, exc):
        message = f"Error while processing {item}: {exc}"
        if self.batchErrorPolicy == "prompt" and not askYesNo(
            f"{message}\nSkip this item and continue? [y/n] "
        ):
            raise exc
        self.errors.append(
            {"item": str(item), "error": str(exc), "exception": exc}
        )
        warnings.warn(f"{message}; skipping it.", RuntimeWarning)

    def _measureSingle(self, runDir, arm, EL, runNumber, rot):
        return measureRunMotorMap(
            runDir,
            arm,
            EL,
            runNumber=runNumber,
            rot=rot,
            xml=self.xml,
            thetaGeometryHdf5=self.thetaGeometryHdf5,
            phiGeometryHdf5=self.phiGeometryHdf5,
            motorMapHdf5=self.motorMapHdf5,
            dataRoot=self.dataRoot,
            binSizeDegrees=self.binSizeDegrees,
            regions=self.regions,
            delta=self.delta,
            steps=self.steps,
            update=self.update,
        )

    def _measureRunDirFile(self, filename):
        table = np.loadtxt(filename, dtype=str, ndmin=2)
        if table.shape[1] != 5:
            raise ValueError(
                f"{filename} must contain: "
                "runNumber, runDir, arm, EL, and rot."
            )

        count = len(table)
        runNumbers = np.full(count, None, dtype=object)
        runDirs = np.full(count, None, dtype=object)
        arms = np.full(count, None, dtype=object)
        elevations = np.full(count, np.nan, dtype=float)
        rotations = np.full(count, np.nan, dtype=float)
        for index, row in enumerate(table):
            try:
                runNumber = int(row[0])
                runDir = str(row[1])
                runDir2Date(runDir)
                arm = validateArm(row[2])
                EL = normalizeEL(row[3])
                rot = normalizeEL(row[4])
                self._measureSingle(runDir, arm, EL, runNumber, rot)
                runNumbers[index] = runNumber
                runDirs[index] = runDir
                arms[index] = arm
                elevations[index] = EL
                rotations[index] = rot
            except Exception as exc:
                self._handleBatchError(
                    f"row {index + 1} ({' '.join(row)})", exc
                )

        self.runNumbers = runNumbers
        self.runDirs = runDirs
        self.arms = arms
        self.ELs = elevations
        self.rots = rotations
        return runDirs

    def determineMotorMap(self):
        if self.runDir is not None:
            runDirInput = str(self.runDir)
            try:
                runDir2Date(runDirInput)
                single = True
            except ValueError:
                single = False

            if single:
                self._validateSingleInput()
                self.runDir = self._measureSingle(
                    runDirInput,
                    self.arm,
                    self.EL,
                    self.runNumber,
                    self.rot,
                )
                return self.runDir
            return self._measureRunDirFile(runDirInput)

        if self.visit is None:
            raise ValueError("Neither visit nor runDir is defined.")

        if isinstance(self.visit, (int, np.integer)):
            self._validateSingleInput()
            self.runDir = visit2runDir(
                int(self.visit),
                self.arm,
                self.EL,
                visitRunDirList=self.visitRunDirList,
                update=self.update,
            )
            self._measureSingle(
                self.runDir,
                self.arm,
                self.EL,
                self.runNumber,
                self.rot,
            )
            return self.runDir

        if isinstance(self.visit, (str, pathlib.Path)):
            (
                self.runNumbers,
                self.runDirs,
                self.visits,
                self.arms,
                self.ELs,
                self.rots,
            ) = visitList2runDirList(
                str(self.visit),
                visitRunDirList=self.visitRunDirList,
                update=self.update,
                errorHandler=self._handleBatchError,
            )
            for runNumber, runDir, arm, EL, rot in zip(
                self.runNumbers,
                self.runDirs,
                self.arms,
                self.ELs,
                self.rots,
            ):
                if runDir is None:
                    continue
                try:
                    self._measureSingle(runDir, arm, EL, runNumber, rot)
                except Exception as exc:
                    self._handleBatchError(
                        f"runDir {runDir} (arm={arm}, EL={EL})", exc
                    )
            return self.runDirs

        raise TypeError("visit must be an integer or a filename.")

    def combineMotorMaps(self):
        return mergeMotorMapsByRunELRotArm(
            self.motorMapHdf5,
            self.combinedMotorMapHdf5,
            overwrite=self.overwriteCombined,
        )

    def run(self, combine=True):
        result = self.determineMotorMap()
        if combine:
            self.combineMotorMaps()
        return result


motorMapMeasurement = MotorMapMeasurement


def parseVisit(value):
    try:
        return int(value)
    except ValueError:
        return value


def buildParser():
    parser = argparse.ArgumentParser(
        description=(
            "Measure run-level PFS cobra motor maps and median-combine "
            "them by run number, elevation, rotation, and arm."
        )
    )
    inputs = parser.add_mutually_exclusive_group(required=True)
    inputs.add_argument("-visit", "--visit", type=parseVisit)
    inputs.add_argument("-runDir", "--runDir")
    parser.add_argument("-arm", "--arm", choices=("theta", "phi"))
    parser.add_argument("-EL", "--EL", type=float)
    parser.add_argument("-runNumber", "--runNumber", "-run", type=int)
    parser.add_argument("-rot", "--rot", type=float)
    parser.add_argument("--xml", default=DEFAULT_XML)
    parser.add_argument(
        "--visit-runDir-list", default="visit_runDir.csv"
    )
    parser.add_argument(
        "--theta-geometry-hdf5",
        default="geometry/thetaGeometry.hdf5",
    )
    parser.add_argument(
        "--phi-geometry-hdf5",
        default="geometry/phiGeometry.hdf5",
    )
    parser.add_argument(
        "--motor-map-hdf5", default="motorMaps/motorMaps.hdf5"
    )
    parser.add_argument(
        "--combined-motor-map-hdf5",
        default="motorMaps/combinedMotorMaps.hdf5",
    )
    parser.add_argument("--data-root", default="/data/MCS")
    parser.add_argument("--bin-size-degrees", type=float, default=3.6)
    parser.add_argument("--regions", type=int, default=112)
    parser.add_argument("--delta", type=float, default=0.1)
    parser.add_argument("--steps", type=int)
    parser.add_argument("--no-update", action="store_true")
    parser.add_argument("--no-combine", action="store_true")
    parser.add_argument("--no-overwrite-combined", action="store_true")
    return parser


def main():
    parser = buildParser()
    args = parser.parse_args()

    scalarInput = isinstance(args.visit, int)
    if args.runDir is not None:
        try:
            runDir2Date(args.runDir)
            scalarInput = True
        except ValueError:
            scalarInput = False
    if scalarInput and any(
        value is None
        for value in (args.arm, args.EL, args.runNumber, args.rot)
    ):
        parser.error(
            "-arm, -EL, -runNumber, and -rot are required "
            "for a single input."
        )

    measurement = MotorMapMeasurement(
        visit=args.visit,
        runDir=args.runDir,
        arm=args.arm,
        EL=args.EL,
        runNumber=args.runNumber,
        rot=args.rot,
        xml=args.xml,
        visitRunDirList=args.visit_runDir_list,
        thetaGeometryHdf5=args.theta_geometry_hdf5,
        phiGeometryHdf5=args.phi_geometry_hdf5,
        motorMapHdf5=args.motor_map_hdf5,
        combinedMotorMapHdf5=args.combined_motor_map_hdf5,
        dataRoot=args.data_root,
        binSizeDegrees=args.bin_size_degrees,
        regions=args.regions,
        delta=args.delta,
        steps=args.steps,
        update=not args.no_update,
        overwriteCombined=not args.no_overwrite_combined,
        batchErrorPolicy="prompt",
    )
    measurement.run(combine=not args.no_combine)


if __name__ == "__main__":
    main()
