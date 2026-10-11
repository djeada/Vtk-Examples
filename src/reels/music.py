"""Generative, tempo-locked soundtrack for science reels.

The track is synthesized from scratch with NumPy: a chord progression drives a
detuned pad, an FM pluck arpeggio, a saturated bass, a lead motif, and a small
drum kit, mixed through a sidechain, ping-pong delay, and convolution reverb.
The tempo is chosen so each story chapter spans a whole number of bars and the
arrangement builds chapter by chapter, so musical changes land on caption
changes. The same seed always produces the same audio.
"""

from __future__ import annotations

from functools import lru_cache
import math

import numpy as np


SAMPLE_RATE = 48_000
_TABLE_SIZE = 4096

_MODES = {
    "minor": (0, 2, 3, 5, 7, 8, 10),
    "dorian": (0, 2, 3, 5, 7, 9, 10),
    "major": (0, 2, 4, 5, 7, 9, 11),
}
# Scale-degree chord roots, one chord per bar.
_PROGRESSIONS = {
    "minor": ((0, 5, 2, 6), (0, 3, 5, 4), (5, 6, 0, 0), (0, 6, 5, 6)),
    "dorian": ((0, 3, 0, 6), (0, 6, 3, 0)),
    "major": ((0, 4, 5, 3), (3, 4, 5, 0), (5, 3, 0, 4)),
}
# Six chapters: intro, groove in, drop, peak, breakdown, final chorus.
_ARRANGEMENT = (
    dict(pad=1.0, arp=0.55, bass=0.0, kick=0, clap=0, hats=0.0, lead=0.0, bright=(0.0, 0.4), impact=False),
    dict(pad=0.9, arp=0.85, bass=0.8, kick=0, clap=0, hats=0.55, lead=0.0, bright=(0.4, 0.6), impact=False),
    dict(pad=0.8, arp=1.0, bass=1.0, kick=1, clap=1, hats=0.8, lead=0.0, bright=(0.6, 0.8), impact=True),
    dict(pad=0.8, arp=1.0, bass=1.0, kick=1, clap=1, hats=1.0, lead=1.0, bright=(0.85, 1.0), impact=True),
    dict(pad=1.0, arp=0.75, bass=0.5, kick=0, clap=0, hats=0.35, lead=0.5, bright=(0.3, 0.55), impact=False),
    dict(pad=0.8, arp=1.0, bass=1.0, kick=1, clap=1, hats=1.0, lead=1.0, bright=(0.8, 1.0), impact=True),
)
_ARP_PATTERNS = (
    (0, 1, 2, 3, 4, 3, 2, 1, 0, 1, 2, 3, 4, 3, 2, 1),
    (0, 2, 1, 3, 2, 4, 3, 5, 0, 2, 1, 3, 2, 4, 3, 5),
    (0, 1, 2, 4, 2, 1, 3, 1, 0, 1, 2, 4, 5, 4, 2, 1),
    (4, 2, 0, 2, 3, 1, 0, 1, 4, 2, 0, 2, 5, 3, 1, 3),
)
# Sixteenth-note onsets per bar.
_BASS_PATTERNS = (
    (0, 3, 6, 8, 11, 14),
    (0, 6, 8, 10, 14),
    (0, 2, 6, 8, 10, 12, 14),
    (0, 7, 8, 11, 14),
)
_LEAD_RHYTHMS = (
    (0, 3, 6, 8, 12, 16, 19, 22, 24),
    (0, 2, 4, 7, 10, 16, 20, 22, 26),
    (0, 4, 6, 10, 12, 16, 18, 24, 28),
)


def _hz(midi):
    return 440.0 * 2.0 ** ((np.asarray(midi, dtype=float) - 69.0) / 12.0)


@lru_cache(maxsize=16)
def _table(kind, max_hz):
    """One band-limited cycle; harmonics stop below Nyquist at ``max_hz``."""
    count = max(1, int(19_000 // max_hz))
    harmonics = np.arange(1, count + 1)
    if kind == "square":
        amplitudes = np.where(harmonics % 2 == 1, 1.0 / harmonics, 0.0)
    else:
        amplitudes = 1.0 / harmonics
    phase = np.arange(_TABLE_SIZE) / _TABLE_SIZE
    wave = amplitudes @ np.sin(2 * np.pi * np.outer(harmonics, phase))
    return wave / np.abs(wave).max()


def _osc(table, frequency, count, phase=0.0):
    """Wavetable oscillator; ``frequency`` is a scalar or per-sample array."""
    if np.ndim(frequency):
        cycles = phase + np.cumsum(frequency) / SAMPLE_RATE
    else:
        cycles = phase + frequency * np.arange(count) / SAMPLE_RATE
    position = (cycles % 1.0) * _TABLE_SIZE
    index = position.astype(np.int64)
    fraction = position - index
    return table[index] * (1 - fraction) + table[(index + 1) % _TABLE_SIZE] * fraction


def _shape(signal, low=None, high=None, order=2):
    """Zero-phase Butterworth-magnitude band shaping in the frequency domain."""
    count = signal.shape[0]
    pad = min(count, SAMPLE_RATE // 2)
    size = count + pad
    spectrum = np.fft.rfft(signal, size, axis=0)
    frequency = np.fft.rfftfreq(size, 1 / SAMPLE_RATE)
    gain = np.ones_like(frequency)
    if high is not None:
        gain /= np.sqrt(1 + (frequency / high) ** (2 * order))
    if low is not None:
        gain /= np.sqrt(1 + (low / np.maximum(frequency, 1e-6)) ** (2 * order))
    if signal.ndim == 2:
        gain = gain[:, None]
    return np.fft.irfft(spectrum * gain, size, axis=0)[:count]


def _convolve(signal, response):
    size = 1 << math.ceil(math.log2(signal.shape[0] + response.shape[0] - 1))
    spectrum = np.fft.rfft(signal, size, axis=0) * np.fft.rfft(response, size, axis=0)
    return np.fft.irfft(spectrum, size, axis=0)[: signal.shape[0]]


def _add(buffer, start, signal, gain=1.0, pan=0.0):
    begin = int(round(start * SAMPLE_RATE))
    if begin >= len(buffer) or gain == 0:
        return
    end = min(len(buffer), begin + len(signal))
    piece = signal[: end - begin] * gain
    if piece.ndim == 1:
        angle = (pan + 1) * np.pi / 4
        buffer[begin:end, 0] += piece * math.cos(angle)
        buffer[begin:end, 1] += piece * math.sin(angle)
    else:
        buffer[begin:end] += piece


def _time(count):
    return np.arange(count) / SAMPLE_RATE


def _pad_note(hz, seconds, attack=0.45, release=0.9):
    count = int((seconds + release) * SAMPLE_RATE)
    table = _table("saw", float(hz) * 1.01)
    voices = np.zeros((count, 2))
    for cents, pan in ((-9, -0.7), (0, 0.0), (8, 0.7)):
        tone = _osc(table, hz * 2 ** (cents / 1200), count, phase=(cents + 9) / 37)
        angle = (pan + 1) * np.pi / 4
        voices[:, 0] += tone * math.cos(angle)
        voices[:, 1] += tone * math.sin(angle)
    t = _time(count)
    envelope = np.minimum(1, t / attack) ** 2
    tail = int(seconds * SAMPLE_RATE)
    envelope[tail:] *= np.cos(np.linspace(0, np.pi / 2, count - tail)) ** 2
    return voices * envelope[:, None] / 3


def _pluck(hz, brightness, seconds=0.45):
    t = _time(int(seconds * SAMPLE_RATE))
    index = (0.5 + 2.4 * brightness) * np.exp(-t * 11)
    tone = np.sin(2 * np.pi * hz * t + index * np.sin(2 * np.pi * hz * 2 * t))
    tone += 0.25 * np.sin(2 * np.pi * hz * 3 * t) * np.exp(-t * 30) * brightness
    return tone * np.exp(-t * 7.5) * np.minimum(1, t / 0.003)


def _bass(hz, seconds):
    count = int((seconds + 0.04) * SAMPLE_RATE)
    t = _time(count)
    tone = (
        np.sin(2 * np.pi * hz * t)
        + 0.35 * np.sin(4 * np.pi * hz * t)
        + 0.12 * np.sin(6 * np.pi * hz * t)
    )
    tone = np.tanh(1.8 * tone) / np.tanh(1.8)
    envelope = np.minimum(1, t / 0.006) * (0.72 + 0.28 * np.exp(-t * 9))
    release = int(0.04 * SAMPLE_RATE)
    envelope[-release:] *= np.linspace(1, 0, release)
    return tone * envelope


def _lead(hz, seconds):
    count = int((seconds + 0.12) * SAMPLE_RATE)
    t = _time(count)
    vibrato = 1 + (2 ** (0.18 / 12) - 1) * np.sin(2 * np.pi * 5.2 * t) * np.clip(
        (t - 0.18) / 0.25, 0, 1
    )
    tone = 0.6 * _osc(_table("square", float(hz) * 1.05), hz * vibrato, count)
    tone += 0.4 * _osc(_table("saw", float(hz) * 1.05), hz * vibrato * 1.003, count)
    envelope = np.minimum(1, t / 0.012) * (0.7 + 0.3 * np.exp(-t * 6))
    release = int(0.12 * SAMPLE_RATE)
    envelope[-release:] *= np.linspace(1, 0, release) ** 2
    return tone * envelope


class _Kit:
    """Fixed drum samples; noise comes from the seeded generator."""

    def __init__(self, rng):
        t = _time(int(0.5 * SAMPLE_RATE))
        sweep = 46 + 120 * np.exp(-t * 30)
        body = np.sin(2 * np.pi * np.cumsum(sweep) / SAMPLE_RATE) * np.exp(-t * 7)
        click = rng.standard_normal(len(t)) * np.exp(-t * 900) * 0.35
        self.kick = np.tanh(1.5 * (body + click)) * 0.95

        t = _time(int(0.4 * SAMPLE_RATE))
        noise = _shape(rng.standard_normal(len(t)), low=900, high=6500)
        bursts = sum(
            np.where(t >= start, np.exp(-(t - start) * 240), 0)
            for start in (0, 0.011, 0.023)
        )
        clap = noise * (bursts + 0.7 * np.exp(-t * 13))
        clap += 0.5 * np.sin(2 * np.pi * 185 * t) * np.exp(-t * 28)
        self.clap = clap / np.abs(clap).max()

        hat_noise = _shape(rng.standard_normal(int(0.3 * SAMPLE_RATE)), low=7000, order=3)
        t = _time(len(hat_noise))
        self.hat = hat_noise * np.exp(-t / 0.035)
        self.hat /= np.abs(self.hat).max()
        self.open_hat = hat_noise * np.exp(-t / 0.16)
        self.open_hat /= np.abs(self.open_hat).max()

        t = _time(int(3.0 * SAMPLE_RATE))
        crash = _shape(rng.standard_normal((len(t), 2)), low=3500) * np.exp(-t / 0.9)[:, None]
        boom = np.sin(2 * np.pi * np.cumsum(38 + 40 * np.exp(-t * 8)) / SAMPLE_RATE)
        crash += 1.4 * (boom * np.exp(-t * 2.8))[:, None]
        self.crash = crash / np.abs(crash).max()

        self.noise = rng.standard_normal((int(8 * SAMPLE_RATE), 2))

    def riser(self, seconds):
        count = int(seconds * SAMPLE_RATE)
        if count < 16:
            return np.zeros((0, 2))
        t = np.linspace(0, 1, count)
        noise = self.noise[np.arange(count) % len(self.noise)]
        bright = _shape(noise, low=2500, high=11000)
        dark = _shape(noise, low=400, high=2500)
        swell = bright * t[:, None] ** 1.5 + dark * (1 - t[:, None]) * 0.5
        tone = np.sin(2 * np.pi * np.cumsum(220 * 4 ** t) / SAMPLE_RATE) * 0.25
        signal = (swell + tone[:, None]) * (t**2.2)[:, None]
        fade = min(count, int(0.01 * SAMPLE_RATE))
        signal[-fade:] *= np.linspace(1, 0, fade)[:, None]
        return signal / max(1e-9, np.abs(signal).max())


def _reverb(signal, rng, seconds=2.4, predelay=0.024):
    t = _time(int(seconds * SAMPLE_RATE))
    response = rng.standard_normal((len(t), 2)) * np.exp(-6.9 * t / seconds)[:, None]
    # High frequencies decay faster: blend a dark full tail with a short bright one.
    bright = np.exp(-6.9 * t / (seconds * 0.35))[:, None]
    response = _shape(response, high=3200) * (1 - bright) + response * bright * 0.6
    response = np.vstack((np.zeros((int(predelay * SAMPLE_RATE), 2)), response))
    response /= np.sqrt((response**2).sum(axis=0))
    return _convolve(signal, response)


def _pingpong(signal, delay, feedback=0.4, taps=6):
    offset = int(delay * SAMPLE_RATE)
    mono = signal.mean(axis=1)
    output = np.zeros_like(signal)
    for tap in range(1, taps + 1):
        shift = tap * offset
        if shift >= len(signal):
            break
        output[shift:, tap % 2] += mono[:-shift] * feedback**tap
    return _shape(output, low=250, high=4500)


def _duck(count, times, depth=0.6, release=0.16):
    envelope = np.ones(count)
    length = int(release * 4 * SAMPLE_RATE)
    curve = 1 - depth * np.exp(-np.arange(length) / (release * SAMPLE_RATE))
    attack = int(0.004 * SAMPLE_RATE)
    curve[:attack] = np.linspace(1, curve[attack], attack)
    for time in times:
        begin = int(round(time * SAMPLE_RATE))
        piece = envelope[begin : begin + length]
        np.minimum(piece, curve[: len(piece)], out=piece)
    return envelope


def _tempo(duration, sections, target):
    """Pick a 4/4 tempo near ``target`` whose bars tile each chapter exactly."""
    chapter = duration / sections
    candidates = [
        240 * bars / chapter
        for bars in range(1, 129)
        if 84 <= 240 * bars / chapter <= 128
    ]
    return min(candidates, key=lambda bpm: abs(bpm - target)) if candidates else target


def render_track(duration, seed=0, sections=6):
    """Return a float stereo track of shape ``(samples, 2)`` peaking near -1 dBFS."""
    rng = np.random.default_rng(seed)
    count = max(1, round(duration * SAMPLE_RATE))
    sections = max(1, int(sections))
    chapter = duration / sections
    mode = ("minor", "dorian", "major")[int(rng.integers(3))]
    scale = _MODES[mode]
    progressions = _PROGRESSIONS[mode]
    progression = progressions[int(rng.integers(len(progressions)))]
    tonic = 45 + int(rng.integers(8))  # A2 .. E3
    bpm = _tempo(duration, sections, float(rng.uniform(98, 118)))
    beat = 60 / bpm
    bar = 4 * beat
    step = beat / 4
    swing = 0.12 * step
    arp_pattern = _ARP_PATTERNS[int(rng.integers(len(_ARP_PATTERNS)))]
    bass_pattern = _BASS_PATTERNS[int(rng.integers(len(_BASS_PATTERNS)))]
    lead_rhythm = _LEAD_RHYTHMS[int(rng.integers(len(_LEAD_RHYTHMS)))]
    lead_tones = rng.integers(0, 5, size=len(lead_rhythm))
    lead_tones[0] = 2
    kit = _Kit(rng)

    def degree(index, octave=0):
        return tonic + 12 * (octave + index // 7) + scale[index % 7]

    def chord(root):
        return [degree(root + offset) for offset in (0, 2, 4, 6)]

    def arrangement(time):
        section = min(sections - 1, int(time / chapter + 1e-9)) if chapter > 0 else 0
        mapped = round(section * (len(_ARRANGEMENT) - 1) / max(1, sections - 1))
        return section, _ARRANGEMENT[mapped]

    pad = np.zeros((count, 2))
    arp = np.zeros((count, 2))
    bass = np.zeros((count, 2))
    lead = np.zeros((count, 2))
    drums = np.zeros((count, 2))
    fx = np.zeros((count, 2))
    clap_send = np.zeros((count, 2))
    kicks = []
    bars = math.ceil(duration / bar - 1e-9)
    final_bar = bars - 1 if bars >= 8 else None

    for index in range(bars):
        start = index * bar
        _, part = arrangement(start + 1e-6)
        root = progression[index % len(progression)]
        if index == final_bar:
            root = 0
        tones = chord(root)
        # Voice the pad around middle C so chord changes move smoothly.
        voiced = sorted(60 + ((tone - 60 + 6) % 12) - 6 for tone in tones)
        for tone in voiced:
            seconds = bar * (2 if index == final_bar else 1)
            _add(pad, start, _pad_note(_hz(tone), seconds), part["pad"] * 0.55)

        if index == final_bar:
            _add(drums, start, kit.kick, 0.9)
            _add(fx, start, kit.crash, 0.35)
            _add(bass, start, _bass(_hz(degree(root, -1)), bar * 1.5), 0.8)
            kicks.append(start)
            continue

        arp_notes = [tone + 12 * octave for octave in (1, 2) for tone in tones[:3]]
        for position in range(16):
            time = start + position * step + (swing if position % 2 else 0)
            accent = 1.0 if position % 4 == 0 else 0.72 if position % 2 == 0 else 0.55
            hz = _hz(arp_notes[arp_pattern[position] % len(arp_notes)])
            pan = 0.35 * math.sin(position * math.pi / 4)
            _add(arp, time, _pluck(hz, accent), part["arp"] * accent * 0.32, pan)

        if part["bass"]:
            onsets = bass_pattern if part["kick"] else (0, 8)
            for number, onset in enumerate(onsets):
                following = onsets[number + 1] if number + 1 < len(onsets) else 16
                octave = 0 if onset % 8 else -1
                hz = _hz(degree(root, octave - 1))
                _add(bass, start + onset * step, _bass(hz, (following - onset) * step * 0.92), part["bass"] * 0.55)

        if part["kick"]:
            for quarter in range(4):
                time = start + quarter * beat
                _add(drums, time, kit.kick, 0.9)
                kicks.append(time)
        if part["clap"]:
            for quarter in (1, 3):
                _add(drums, start + quarter * beat, kit.clap, 0.42, 0.05)
                _add(clap_send, start + quarter * beat, kit.clap, 0.42)
        if part["hats"]:
            for position in range(16):
                time = start + position * step + (swing if position % 2 else 0)
                if position % 4 == 2 and part["kick"]:
                    _add(drums, time, kit.open_hat, 0.16 * part["hats"], 0.25)
                elif position % 2 or part["hats"] >= 0.8:
                    velocity = 0.75 if position % 2 else 0.45
                    _add(drums, time, kit.hat, 0.17 * part["hats"] * velocity, 0.3)

        if part["lead"]:
            phrase = index % 2
            for number, onset in enumerate(lead_rhythm):
                if onset // 16 != phrase:
                    continue
                following = lead_rhythm[number + 1] if number + 1 < len(lead_rhythm) else 32
                pitch = tones[lead_tones[number] % 4] + 24 + (12 if lead_tones[number] == 4 else 0)
                seconds = min(following - onset, 6) * step * 0.9
                _add(lead, start + (onset % 16) * step, _lead(_hz(pitch), seconds), part["lead"] * 0.2)

    for section in range(1, sections):
        _, part = arrangement(section * chapter + 1e-6)
        if part["impact"]:
            length = min(bar, chapter)
            _add(fx, section * chapter - length, kit.riser(length), 0.22)
            _add(fx, section * chapter, kit.crash, 0.3)

    t = np.arange(count) / SAMPLE_RATE
    times, values = [], []
    for section in range(sections):
        _, part = arrangement(section * chapter + 1e-6)
        times += [section * chapter, (section + 1) * chapter - 1e-3]
        values += list(part["bright"])
    brightness = np.interp(t, times, values)[:, None]
    pad = _shape(pad, high=650) * (1 - brightness) + _shape(pad, high=4200) * brightness
    arp = _shape(arp, high=1400) * (1 - brightness) + _shape(arp, high=8000) * brightness
    lead = _shape(lead, low=180, high=5200)
    bass = _shape(bass, high=1400)

    duck = _duck(count, kicks)[:, None]
    dry = pad * duck + arp * (0.55 + 0.45 * duck) + bass * duck + lead + drums + fx
    sends = 0.5 * pad + 0.45 * arp + 0.4 * lead + 0.25 * clap_send + 0.15 * fx
    echoes = _pingpong(arp * 0.5 + lead, 0.75 * beat)
    wet = _reverb(sends + 0.4 * echoes, rng) * duck
    mix = dry + 0.32 * echoes + 0.5 * wet
    mix = _shape(mix, low=28)

    fade = min(2.5, duration * 0.2)
    if fade > 0:
        length = min(count, int(fade * SAMPLE_RATE))
        mix[-length:] *= (np.cos(np.linspace(0, np.pi / 2, length)) ** 2)[:, None]
    lead_in = min(count, int(0.02 * SAMPLE_RATE))
    mix[:lead_in] *= np.linspace(0, 1, lead_in)[:, None]

    peak = np.abs(mix).max()
    if peak > 0:
        mix /= peak
    # Gentle glue saturation, then leave 1 dB of headroom for the AAC encoder.
    mix = np.tanh(1.4 * mix) / np.tanh(1.4)
    return mix * 10 ** (-1 / 20)
