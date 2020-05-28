def test(**kwargs):
    for a in kwargs:
        if type(kwargs[a]) is list:
            for b in kwargs[a]:
                print(b)
        else :
            print(kwargs[a])

test(a = [1, 2, 3, 4, 5, 6], b = 22, c = 33)

symbol = 'aaaaa'

eval(symbol) = 42

print(aaaaa)