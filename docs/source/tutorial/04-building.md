# building docs locally

## venv
```sh
python3 -m venv .venv
```
enter the venv
````{tab} *nux
```sh
source .venv/bin/activate
```
````
````{tab} windows
```sh
.venv/Scripts/activate
```
````
install required packages
```sh
pip install -r requirements
```

## manual build
````{tab} *nux
```sh
make html
```
````
````{tab} windows
```{caution}
not tested
```
```sh
make.bat html
```
````

## auto rebuilding
```sh
sphinx-autobuild source/ build/
```