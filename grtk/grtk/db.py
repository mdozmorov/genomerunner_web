import psycopg2

def connect():
    return psycopg2.connect(database="gfa", host="wrendb", user="gilesc")
