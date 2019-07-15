import click

@click.command()
@click.option('--hello')
@click.argument('name')
def test(hello,name):
    print(f"hello: {hello}")
    print(f"name: {name}")

if __name__ == "__main__":
    test()
