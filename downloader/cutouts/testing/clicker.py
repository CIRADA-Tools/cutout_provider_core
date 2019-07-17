import click

@click.command()
@click.option('--hello')
@click.argument('name')
#@click.option('--bart','simpsons',is_flag=True,default=None,cls='MutuallyExclusiveOption')
#@click.option('--lisa','simpsons',is_flag=True,default=None,cls='MutuallyExclusiveOption')
def test(hello,name,simpsons):
    print(f"hello: {hello}")
    print(f"name: {name}")
    print(f"foo: {simpsons}")

if __name__ == "__main__":
    test()
