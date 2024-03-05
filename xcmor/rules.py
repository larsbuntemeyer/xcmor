from .log import get_logger
from .mapping import dtype_map

logger = get_logger(__name__)


class rules:
    # attributes are dropped after interpretation
    drop = True

    @classmethod
    def type(cls, obj, time=False):
        """apply dtype rule to dataarray"""
        if time is False:
            axis = obj.attrs.get("axis")
            # don't convert time variables.
            if axis and axis == "T":
                return obj
        dtype = obj.attrs["type"]
        if obj.dtype != dtype_map[dtype]:
            logger.warning(
                f"converting {obj.name or 'data'} from {obj.dtype} to {dtype_map[dtype]}"
            )
            obj = obj.astype(dtype_map[dtype])
        if cls.drop is True:
            del obj.attrs["type"]
        return obj

    @classmethod
    def out_name(cls, obj):
        """apply renaming rule"""
        out_name = obj.attrs["out_name"]
        obj.name = out_name
        # obj.rename_vars({v: out_name})
        if cls.drop is True:
            del obj.attrs["out_name"]
        return obj

    @classmethod
    def valid_min(cls, obj):
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
        return obj

    @classmethod
    def valid_max(cls, obj):
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
        return obj

    @classmethod
    def standard_name(cls, obj):
        from cf_xarray.utils import parse_cf_standard_name_table

        sname = obj.attrs["standard_name"]
        info, table, aliases = parse_cf_standard_name_table()
        if sname not in table.keys():
            raise Exception(
                f"'{obj.name+'': ' or ''}{sname} is not a valid standard name."
            )
        return obj

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
