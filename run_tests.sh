#!/bin/bash
export CONFIG_OBJ='config.TestingConfig'
nosetests-2.7 -v
unset CONFIG_OBJ
