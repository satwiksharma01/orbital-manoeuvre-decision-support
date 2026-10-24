"""Record the web planner's animations into the README demo GIF.

    python tools/record_demo_gif.py [page.html] [out.gif] [--still out.png]

Defaults: web/index.html -> assets/demo.gif. Needs Microsoft Edge, `pip install selenium pillow`.
Headless Edge at 1280 x 1000 in dark mode clicks four example missions. The page's clock runs
4x slow while recording, so a screenshot per frame keeps up with the animations; frames are
played back at real speed. --still also saves a PNG of the LEO -> Molniya state."""
import base64
import io
import sys
import time
from pathlib import Path

from PIL import Image
from selenium import webdriver
from selenium.webdriver.common.by import By

ROOT = Path(__file__).resolve().parent.parent
still = sys.argv[sys.argv.index("--still") + 1] if "--still" in sys.argv else None
pos = [a for a in sys.argv[1:] if a != "--still" and a != still]
page = Path(pos[0] if pos else ROOT / "web" / "index.html").resolve().as_uri()
out = pos[1] if len(pos) > 1 else str(ROOT / "assets" / "demo.gif")

opts = webdriver.EdgeOptions()
opts.add_argument("--headless=new")
opts.add_argument("--window-size=1280,1000")
opts.add_argument("--hide-scrollbars")
opts.add_argument("--force-dark-mode")
d = webdriver.Edge(options=opts)
d.set_window_size(1280, 1000)
SLOW = 4   # page clock runs 4x slower while recording; frames play back at real speed
d.execute_cdp_cmd("Page.addScriptToEvaluateOnNewDocument", {"source": f"""
  const pn = performance.now.bind(performance), raf = window.requestAnimationFrame.bind(window);
  performance.now = () => pn() / {SLOW};
  window.requestAnimationFrame = cb => raf(t => cb(t / {SLOW}));
"""})
d.get(page)
time.sleep(1.5)
d.execute_script("window.scrollTo(0, 250)")   # Orbits card heading just below the toolbar
time.sleep(0.5)

frames, durations = [], []


def capture(seconds):
    end = time.time() + seconds
    while time.time() < end:
        t0 = time.time()
        shot = d.execute_cdp_cmd("Page.captureScreenshot", {"format": "jpeg", "quality": 85})
        img = Image.open(io.BytesIO(base64.b64decode(shot["data"]))).convert("RGB")
        frames.append(img.resize((960, int(img.height * 960 / img.width)), Image.LANCZOS))
        durations.append((time.time() - t0) / SLOW)


def click(label):
    for b in d.find_elements(By.CSS_SELECTOR, ".chip"):
        if b.text.strip() == label:
            b.click()
            return
    raise SystemExit(f"chip not found: {label}")


capture(1.2 * SLOW)
for label in ["36-min rapid response", "GTO → GEO", "LEO → Molniya", "LEO → GEO"]:
    click(label)
    capture(2.2 * SLOW)
if still:
    click("LEO → Molniya")
    time.sleep(3 * SLOW)
    d.execute_cdp_cmd("Page.captureScreenshot", {"format": "png"})   # warm-up
    shot = d.execute_cdp_cmd("Page.captureScreenshot", {"format": "png"})
    Image.open(io.BytesIO(base64.b64decode(shot["data"]))).save(still)
d.quit()

# keep real-time playback: merge frames so each shown frame lasts >= 40 ms (GIF players clamp shorter delays)
kept, ms, acc = [], [], 0.0
for f, s in zip(frames, durations):
    if acc == 0.0:
        kept.append(f)
    acc += s
    if acc >= 0.04:
        ms.append(int(round(acc * 1000))); acc = 0.0
if len(ms) < len(kept):
    ms.append(max(40, int(round(acc * 1000))))
frames = kept
ms[-1] = 1500                              # hold the last frame before looping
sample = frames[:: max(1, len(frames) // 12)]       # palette from frames across all missions
SERIES = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#1f9fb4", "#8a63d2",    # light-theme series
          "#3987e5", "#d95926", "#199e70", "#c98500", "#d55181", "#1b93a7", "#9873e3"]    # dark-theme series
sw = 120                                                          # swatches keep rare series colours
strip = Image.new("RGB", (frames[0].width, frames[0].height * len(sample) + sw))
for i, f in enumerate(sample):
    strip.paste(f, (0, i * f.height))
for i, c in enumerate(SERIES):
    strip.paste(c, (i * strip.width // len(SERIES), strip.height - sw, (i + 1) * strip.width // len(SERIES), strip.height))
pal = strip.quantize(colors=256, method=Image.Quantize.MEDIANCUT)
q = [f.quantize(palette=pal, dither=Image.Dither.NONE) for f in frames]
q[0].save(out, save_all=True, append_images=q[1:], duration=ms, loop=0, optimize=True)
print(f"{len(frames)} frames, {sum(ms)/1000:.1f} s, {Path(out).stat().st_size/1e6:.1f} MB")
