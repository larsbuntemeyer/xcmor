from .log import get_logger
from .mapping import dtype_map

logger = get_logger(__name__)


class rules:
    # attributes are dropped after interpretation
    drop = True

    @classmethod
    def type(cls, ds, name):
        obj = ds[name]
        """apply dtype rule to dataarray"""
        axis = obj.attrs.get("axis")
        dtype = obj.attrs["type"]
        # don't convert time variables, just add encoding
        if axis and axis == "T":
            obj.encoding["dtype"] = dtype_map[dtype]
            del obj.attrs["type"]
            return ds
        if obj.dtype != dtype_map[dtype]:
            logger.info(
                f"converting {obj.name or 'data'} from {obj.dtype} to {dtype_map[dtype]}"
            )
            ds[name] = obj.astype(dtype_map[dtype])

        return ds

    @classmethod
    def out_name(cls, ds, name):
        """apply renaming rule"""
        obj = ds[name]
        out_name = obj.attrs["out_name"]
        obj.name = out_name
        # obj.rename_vars({v: out_name})
        if cls.drop is True:
            del obj.attrs["out_name"]
        return ds

    @classmethod
    def valid_min(cls, ds, name):
        obj = ds[name]
        valid_min = obj.attrs.get("valid_min")
        if valid_min:
            try:
                assert obj.min() >= float(valid_min)
            except Exception:
                raise Exception(
                    f"{obj.name or 'data'} is violating valid_min: {valid_min}"
                )
        if cls.drop is True:
            del obj.attrs["valid_min"]
        return ds

    @classmethod
    def valid_max(cls, ds, name):
        obj = ds[name]
        valid_max = obj.attrs.get("valid_max")
        if valid_max:
            try:
                assert obj.min() <= float(valid_max)
            except Exception:
                raise Exception(
                    f"{obj.name or 'data'} is violating valid_max: {valid_max}"
                )
        if cls.drop is True:
            del obj.attrs["valid_max"]
        return ds

    @classmethod
    def standard_name(cls, ds, name):
        obj = ds[name]
        from cf_xarray.utils import parse_cf_standard_name_table

        sname = obj.attrs["standard_name"]
        info, table, aliases = parse_cf_standard_name_table()
        if sname not in table.keys():
            raise Exception(
                f"'{obj.name+'': ' or ''}{sname} is not a valid standard name."
            )
        return ds

    @classmethod
    def must_have_bounds(cls, ds, name):
        obj = ds[name]
        bounds_required = obj.attrs.get("must_have_bounds") == "yes"
        if bounds_required:
            logger.debug(f"checking bounds for {name}: {bounds_required}")
        if bounds_required and not ds.cf.bounds.get(name):
            logger.warning(f"{name} must have bounds")
            try:
                logger.info(f"adding bounds for {name}")
                return ds.cf.add_bounds(name)
            except Exception as e:
                logger.error(f"Failed to add bounds to {name}: {e}")
        return ds

    # @classmethod
    # def requested(cls, obj):
    #     requested = obj.attrs["requested"]
    #     if requested:
    #         requested = list(map(float, requested))
    #     if not np.allclose(obj.values, requested):
    #         raise Exception(
    #             f"{obj.name or ' '} has not all inconsistent values: {requested} "
    #         )

    # @classmethod
    # def frequency(cls, obj):
    #     raise NotImplementedError

    # @classmethod
    # def ok_min_mean_abs(cls, obj):
    #     raise NotImplementedError

    # @classmethod
    # def ok_max_mean_abs(cls, obj):
    #     raise NotImplementedError

    # @classmethod
    # def cell_measures(cls, obj):
    #     raise NotImplementedError

    # @classmethod
    # def cell_methods(cls, obj):
    #     raise NotImplementedError
