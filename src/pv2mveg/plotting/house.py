"""House style for generated text: ISO months and dates, ΔNDVI for the baseline-adjusted contrast.

Kept free of plotting imports so that any script can use it without changing matplotlib settings."""
import re

_MONTHS = ['January', 'February', 'March', 'April', 'May', 'June', 'July', 'August', 'September', 'October', 'November', 'December']
_MNUM = {m: i + 1 for i, m in enumerate(_MONTHS)} | {m[:3]: i + 1 for i, m in enumerate(_MONTHS)}
_M = '(' + '|'.join(_MONTHS) + ')'
_MA = '(' + '|'.join(_MONTHS + [m[:3] for m in _MONTHS]) + ')'
_DELTA = [
    (r'[Bb]aseline-adjusted (NDVI|NIRv|MSAVI) contrasts?', r'Δ\1'),
    (r'inter-row [Bb]aseline-adjusted contrasts for NDVI, NIRv and MSAVI', 'inter-row ΔNDVI, ΔNIRv and ΔMSAVI'),
    (r'[Bb]aseline-adjusted (NDVI|NIRv|MSAVI)\b', r'Δ\1'),
    (r'(?<![\w-])[Aa]djusted (NDVI|NIRv|MSAVI)\b', r'Δ\1'),
    (r'the adjusted inter-row NDVI difference Δadj', 'the inter-row ΔNDVI'),
    (r'Adjusted contrast \[joint 95%\]', 'Δ [joint 95%]'),
]


def house(text):
    """ISO months and dates; ΔNDVI for the baseline-adjusted contrast. Idempotent."""
    text = re.sub(r'\b(\d{1,2})–(\d{1,2}) ' + _M + r' (20\d\d)\b',
                  lambda m: f'{m[4]}-{_MNUM[m[3]]:02d}-{int(m[1]):02d} to {m[4]}-{_MNUM[m[3]]:02d}-{int(m[2]):02d}', text)
    text = text.replace('Satellite dates are 17, 13 and 14 August, respectively', 'Satellite dates are 2025-08-17, 2025-08-13 and 2025-08-14, respectively')
    text = re.sub(r'\b(\d{1,2}) ' + _M + r' (20\d\d)\b', lambda m: f'{m[3]}-{_MNUM[m[2]]:02d}-{int(m[1]):02d}', text)
    # every day-month without a year in this paper is a 2025 acquisition
    text = re.sub(r'\b(\d{1,2}) (May|June|July|August)\b(?! 20\d\d)', lambda m: f'2025-{_MNUM[m[2]]:02d}-{int(m[1]):02d}', text)
    text = re.sub(r'\b' + _M + r' − ' + _M + r' (20\d\d)\b', lambda m: f'{m[3]}-{_MNUM[m[1]]:02d} − {m[3]}-{_MNUM[m[2]]:02d}', text)
    text = re.sub(r'\b' + _M + r'-minus-' + _M + r' (20\d\d)\b', lambda m: f'{m[3]}-{_MNUM[m[1]]:02d} minus {m[3]}-{_MNUM[m[2]]:02d}', text)
    text = re.sub(r'\b' + _M + r'-minus-(20\d\d)-(\d\d)\b', lambda m: f'{m[2]}-{_MNUM[m[1]]:02d} minus {m[2]}-{m[3]}', text)
    text = re.sub(r'\b' + _M + r' − (20\d\d)-(\d\d)\b', lambda m: f'{m[2]}-{_MNUM[m[1]]:02d} − {m[2]}-{m[3]}', text)
    text = re.sub(r'\b' + _M + r' (and|to|or) ' + _M + r' (20\d\d)\b', lambda m: f'{m[4]}-{_MNUM[m[1]]:02d} {m[2]} {m[4]}-{_MNUM[m[3]]:02d}', text)
    text = re.sub(r'\b' + _MA + r'\.? (20\d\d)\b', lambda m: f'{m[2]}-{_MNUM[m[1]]:02d}', text)
    text = re.sub(r'(20\d\d-\d\d)-minus-(20\d\d-\d\d)', r'\1 minus \2', text)
    for a, b in _DELTA:
        text = re.sub(a, b, text)
    return text
