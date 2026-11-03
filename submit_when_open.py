#!/usr/bin/env python
"""Verify a CYP challenge submission now, wait, then submit.

Run detached so a closing terminal cannot kill it:

    conda run -n <env-with-gradio_client> python submit_when_open.py      # tmux
    nohup python submit_when_open.py > /dev/null 2>&1 &                   # nohup

Everything checkable is checked BEFORE the wait, while you are still at the
keyboard. The 12 h lockout makes a bad submission cost the window as well as
the attempt. Progress goes to both stdout and a log file, so `tail -f` works
whichever way you launched it.
"""
import os, sys, time, zipfile, socket
from pathlib import Path
from datetime import datetime, timedelta
from zoneinfo import ZoneInfo

from dotenv import load_dotenv, find_dotenv
import huggingface_hub
from huggingface_hub import get_token, whoami
from gradio_client import Client, handle_file     # <- missing from the cell

# ---------------------------------------------------------------- config ---
TZ          = ZoneInfo("America/New_York")   # survives the Nov 2 DST change
SUBMIT_AT   = "2026-11-01 04:10"             # local time; None = submit now
EXPECT_N    = 20                             # files expected in the zip
POLL_MIN    = 15                             # heartbeat interval, minutes
RETRY_MIN   = 15                             # wait between retries on lockout
MAX_RETRIES = 6
DRY_RUN     = False                          # True = preflight + wait, no send
LOG_DIR     = Path(__file__).resolve().parent
ANON        = True
SPACE       = "openadmet/cyp-challenge"
# ---------------------------------------------------------------------------

STAMP = datetime.now(TZ).strftime("%Y%m%d_%H%M%S")
LOG   = LOG_DIR / f"submit_{STAMP}.log"
STATE = LOG_DIR / "submit_status.txt"


def say(msg=""):
    line = f"{datetime.now(TZ):%H:%M:%S}  {msg}" if msg else ""
    print(line, flush=True)
    with LOG.open("a") as fh:
        fh.write(line + "\n")


def state(msg):
    STATE.write_text(
        f"{datetime.now(TZ):%Y-%m-%d %H:%M:%S %Z}\npid {os.getpid()} on "
        f"{socket.gethostname()}\n{msg}\nlog: {LOG}\n")


def die(msg):
    say(f"[abort] {msg}")
    state(f"ABORTED: {msg}")
    sys.exit(1)


say(f"[start] pid {os.getpid()} on {socket.gethostname()}")
say(f"[start] log {LOG}")
load_dotenv(find_dotenv(), override=True)

# -------------------------------------------------------------- preflight ---
KEYS = ("username", "user_alias", "participant_name", "discord_username",
        "email", "affiliation", "model_tag", "zip_path")
PARAMS = {k: os.getenv(k) for k in KEYS}
missing = [k for k, v in PARAMS.items() if not v]
if missing:
    die(f"unset in .env: {missing}")

zp = Path(PARAMS["zip_path"]).expanduser().resolve()
if not zp.is_file():
    die(f"zip not found: {zp}")

with zipfile.ZipFile(zp) as z:
    names = [n for n in z.namelist() if not n.endswith("/")]
    corrupt = z.testzip()
if corrupt:
    die(f"corrupt member in zip: {corrupt}")

stems = sorted({Path(n).stem for n in names})
say(f"[preflight] {zp}")
say(f"[preflight] {zp.stat().st_size/1e6:.2f} MB, modified "
    f"{datetime.fromtimestamp(zp.stat().st_mtime, TZ):%Y-%m-%d %H:%M %Z}")
say(f"[preflight] {len(names)} files, {len(stems)} unique stems")

problems = []
if len(names) != EXPECT_N:
    problems.append(f"{len(names)} files, expected {EXPECT_N}")
if len(stems) != len(names):
    problems.append("duplicate stems inside the zip")
odd = [s for s in stems if not s.upper().startswith("OCNT")]
if odd:
    problems.append(f"names not in OCNT-* form: {odd[:5]}")
with zipfile.ZipFile(zp) as z:
    for n in names:
        head = z.read(n)[:400].decode("utf-8", "replace")
        if "ATOM" not in head and "HETATM" not in head:
            problems.append(f"{n}: no ATOM/HETATM in the first 400 bytes")
            break
        if " H  " in head:
            problems.append(f"{n}: looks like it still contains hydrogens")
            break

say("[preflight] parameters")
for k, v in PARAMS.items():
    if k == "zip_path":
        continue
    say(f"            {k:18s} {v if k != 'email' else v[:3] + '***'}")
say(f"            {'track':18s} Structure Prediction")
say(f"            {'anonymous':18s} {ANON}      open_code  True")
# ---- auth preflight: everything that must be true at submit time ----------
say("[preflight] auth")
say(f"            huggingface_hub   {huggingface_hub.__version__}")
if tuple(int(x) for x in huggingface_hub.__version__.split(".")[:2]) >= (2, 0):
    problems.append(
        f"huggingface_hub {huggingface_hub.__version__} is incompatible with "
        "gradio_client (needs <2.0) -- pip install 'huggingface_hub==1.29.0'")

TOKEN = os.getenv("HF_TOKEN") or get_token()
if not TOKEN:
    from huggingface_hub import constants as _hfc
    problems.append(
        "no HF token: HF_TOKEN is unset and no token file was found at "
        f"{getattr(_hfc, 'HF_TOKEN_PATH', '~/.cache/huggingface/token')}"
        " -- run `hf auth login`, or put HF_TOKEN in .env")
else:
    try:
        say(f"            authenticated as  {whoami(token=TOKEN).get('name','?')}")
    except Exception as ex:
        problems.append(f"token present but whoami() failed: {str(ex)[:200]}")

    # strongest check available: build the client now, while you are watching
    try:
        _c = Client(SPACE, token=TOKEN, verbose=False)
        if not any("submit_predictions" in str(e) for e in _c.view_api(
                return_format="dict")["named_endpoints"]):
            problems.append("/submit_predictions not found on the Space")
        say(f"            space reachable   {SPACE}")
        del _c
    except Exception as ex:
        problems.append(f"cannot reach {SPACE}: {str(ex)[:200]}")
if problems:
    for p in problems:
        say(f"   [!] {p}")
    die("preflight failed; nothing was sent")
say("[preflight] OK -- model_tag above is what appears on the leaderboard")
if DRY_RUN:
    say("[preflight] DRY_RUN is set: will wait, then stop without submitting")

# ------------------------------------------------------------------ wait ---
if SUBMIT_AT:
    target = datetime.strptime(SUBMIT_AT, "%Y-%m-%d %H:%M").replace(tzinfo=TZ)
    now = datetime.now(TZ)
    if target <= now:
        say(f"[wait] {target:%Y-%m-%d %H:%M} is in the past; going now")
    else:
        say(f"[wait] now {now:%Y-%m-%d %H:%M:%S %Z}; submitting at "
            f"{target:%Y-%m-%d %H:%M %Z} "
            f"({(target-now).total_seconds()/3600:.2f} h)")
        state(f"WAITING until {target:%Y-%m-%d %H:%M %Z}")
        while True:
            left = (target - datetime.now(TZ)).total_seconds()
            if left <= 0:
                break
            time.sleep(min(POLL_MIN * 60, left))
            left = max(0.0, (target - datetime.now(TZ)).total_seconds())
            say(f"[wait] {left/3600:5.2f} h to go")
            state(f"WAITING until {target:%Y-%m-%d %H:%M %Z} "
                  f"({left/3600:.2f} h left)")

# ---------------------------------------------------------------- submit ---
if DRY_RUN:
    say("[submit] DRY_RUN -- stopping here, nothing sent")
    state("DRY RUN complete; nothing submitted")
    sys.exit(0)

try:
    client = Client(SPACE, token=TOKEN)
except Exception as ex:
    die(f"could not reach {SPACE} at submit time: {str(ex)[:300]}")
result = None
for attempt in range(1, MAX_RETRIES + 1):
    try:
        say(f"[submit] attempt {attempt}/{MAX_RETRIES}")
        state(f"SUBMITTING, attempt {attempt}/{MAX_RETRIES}")
        result = client.predict(
            username=PARAMS["username"],
            user_alias=PARAMS["user_alias"],
            anon_checkbox=ANON,
            participant_name=PARAMS["participant_name"],
            discord_username=PARAMS["discord_username"],
            email=PARAMS["email"],
            affiliation=PARAMS["affiliation"],
            model_tag=PARAMS["model_tag"],
            paper_checkbox=True,
            proprietary_data_checkbox=False,
            open_code_checkbox=True,
            track_select="Structure Prediction",
            file_input=handle_file(str(zp)),     # the zip we actually verified
            api_name="/submit_predictions",
        )
        say("[submit] ACCEPTED")
        say(str(result))
        break
    except Exception as ex:
        msg = str(ex)
        say(f"[submit] failed: {msg[:300]}")
        locked = any(w in msg.lower() for w in
                     ("lock", "wait", "too soon", "12 h", "12 hour", "rate"))
        if attempt == MAX_RETRIES:
            die("giving up; nothing submitted")
        if not locked:
            die("not a lockout error -- stopping so you can look")
        say(f"[submit] lockout still active; retrying in {RETRY_MIN} min")
        state(f"LOCKED OUT, retry {attempt+1}/{MAX_RETRIES} "
              f"in {RETRY_MIN} min")
        time.sleep(RETRY_MIN * 60)

done = datetime.now(TZ)
nxt = done + timedelta(hours=12)
say(f"[submit] submitted {done:%Y-%m-%d %H:%M %Z}")
say(f"[submit] next window opens {nxt:%Y-%m-%d %H:%M %Z}"
    f'   -> SUBMIT_AT = "{nxt:%Y-%m-%d %H:%M}"')
state(f"SUBMITTED {done:%Y-%m-%d %H:%M %Z}\n"
      f"model_tag {PARAMS['model_tag']}\nzip {zp.name}\n"
      f"next window {nxt:%Y-%m-%d %H:%M %Z}")
with (LOG_DIR / "submission_log.txt").open("a") as fh:
    fh.write(f"{done:%Y-%m-%d %H:%M %Z}\t{PARAMS['model_tag']}\t{zp.name}\n")
