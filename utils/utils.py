import re
import os
import awkward as ak


def atoi(text):
    return int(text) if text.isdigit() else text


def natural_keys(text):
    """
    alist.sort(key=natural_keys) sorts in human order
    http://nedbatchelder.com/blog/200712/human_sorting.html
    (See Toothy's implementation in the comments)
    """
    return [atoi(c) for c in re.split(r"(\d+)", text)]


def get_files(paths, pattern=".root", exclude=None):
    files = []
    for path in paths:
        for it in os.scandir(path):
            exclude_file = True
            if it.name.find(pattern) > -1 and (it.stat().st_size != 0):
                exclude_file = False
                if exclude:
                    if type(exclude) == str:
                        if it.name.find(exclude) > -1:
                            exclude_file = True
                    else:
                        for e in exclude:
                            if it.name.find(e) > -1:
                                exclude_file = True
            if not exclude_file:
                files.append(it.path)
    files.sort(key=natural_keys)
    return files


def remove_none(arr):
    return arr[~ak.is_none(arr.pt, axis=-1)]


def fill_kin_hists(obj, hists, cat=None):
    for h in hists:
        if h not in obj.fields:
            if h == "n":
                o = ak.num(obj)
            else:
                continue
        else:
            o = ak.flatten(obj[h])

        if cat:
            hists[h].fill(o, cat=cat)
        else:
            hists[h].fill(o)
