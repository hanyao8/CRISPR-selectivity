a=[0,1,1,1,4]
"""
def all_perms(elements):

    for perm in all_perms(elements[1:]):
        for i in range(len(elements)):
                # nb elements[0:1] works in both string and list contexts
            yield perm[:i] + elements[0:1] + perm[i:]

"""
def all_perms(elements):
    if len(elements) <=1:
        yield elements
    else:
        for perm in all_perms(elements[1:]):
            for i in range(len(elements)):
                # nb elements[0:1] works in both string and list contexts
                yield perm[:i] + elements[0:1] + perm[i:]
