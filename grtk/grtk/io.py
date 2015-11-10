import pandas

def read_matrix(path):
    M = pandas.read_table(path, sep="\t", index_col=0)
    M.index = map(str, M.index)
    M.columns = map(str, M.columns)
    return M
